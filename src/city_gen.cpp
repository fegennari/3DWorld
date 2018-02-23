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
float const ROAD_HEIGHT             = 0.001;
float const OUTSIDE_TERRAIN_HEIGHT  = 0.0;


extern int rand_gen_index, display_mode;
extern float water_plane_z, shadow_map_pcf_offset, cobj_z_bias;


struct city_params_t {

	unsigned num_cities, num_samples, city_size, city_border, road_border, slope_width;
	float road_width, road_spacing;

	city_params_t() : num_cities(0), num_samples(100), city_size(0), city_border(0), road_border(0), slope_width(0), road_width(0.0), road_spacing(0.0) {}
	bool enabled() const {return (num_cities > 0 && city_size > 0);}
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
		else if (str == "city_size") {
			if (!read_uint(fp, city_size)) {return read_error(str);}
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

float smooth_interp(float a, float b, float mix) {
	mix = mix * mix * (3.0 - 2.0 * mix); // cubic Hermite interoplation (smoothstep)
	return mix*a + (1.0 - mix)*b;
}

class heightmap_query_t {
protected:
	float *heightmap;
	unsigned xsize, ysize;

public:
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
	void flatten_region_to(cube_t const c, unsigned slope_width) {
		flatten_region_to(get_x_pos(c.x1()), get_y_pos(c.y1()), get_x_pos(c.x2()), get_y_pos(c.y2()), slope_width, (c.z1() - ROAD_HEIGHT));
	}
	void flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation) {
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = max((int)y1-(int)slope_width, 0); y < min(y2+slope_width, ysize); ++y) {
			for (unsigned x = max((int)x1-(int)slope_width, 0); x < min(x2+slope_width, xsize); ++x) {
				float &h(get_height(x, y));

				if (slope_width > 0) {
					float const dx(max(0, max(((int)x1 - (int)x), ((int)x - (int)x2 + 1))));
					float const dy(max(0, max(((int)y1 - (int)y), ((int)y - (int)y2 + 1))));
					h = smooth_interp(h, elevation, min(1.0f, sqrt(dx*dx + dy*dy)/slope_width));
				}
				else {h = elevation;}
			}
		}
	}
	float flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, bool stats_only=0) {
		assert(is_valid_region(x1, y1, x2, y2));
		if (x1 == x2 || y1 == y2) return 0.0; // zero area
		float const run_len(dim ? (y2 - y1) : (x2 - x1)), denom(1.0f/max(run_len, 1.0f)), dz(z2 - z1);
		unsigned const border(city_params.road_border); // just fish it out of the global params rather than pass it all the way down
		int const pad(border + 1U); // pad an extra 1 texel to handle roads misaligned with the texture
		unsigned px1(x1), py1(y1), px2(x2), py2(y2);
		float tot_dz(0.0);
		
		if (dim) {
			px1 = max((int)x1-pad, 0);
			px2 = min(x2+pad, xsize);
			py1 = max((int)y1-1, 0); // pad by 1 in road dim as well to blend with edge of city
			py2 = min(y2+1, ysize);
		}
		else     {
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
				float new_h;

				if (border > 0) {
					unsigned const dist(dim ? max(0, max(((int)x1 - (int)x - 1), ((int)x - (int)x2))) : max(0, max(((int)y1 - (int)y - 1), ((int)y - (int)y2))));
					new_h = smooth_interp(h, road_z, float(dist)/float(border));
				}
				else {new_h = road_z;}
				tot_dz += fabs(h - new_h);
				if (!stats_only) {h = new_h;} // apply the height change
			} // for x
		} // for y
		return tot_dz;
	}
};

class city_plot_gen_t : public heightmap_query_t {

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
	bool find_best_city_location(unsigned width, unsigned height, unsigned border, unsigned slope_width, unsigned num_samples, unsigned &x_llc, unsigned &y_llc) {
		assert(num_samples > 0);
		assert((width + 2*border) < xsize && (height + 2*border) < ysize); // otherwise the city can't fit in the map
		unsigned const num_iters(100*num_samples); // upper bound
		unsigned xend(xsize - width - 2*border + 1), yend(ysize - width - 2*border + 1); // max rect LLC, inclusive
		unsigned num_cands(0);
		float best_diff(0.0);

		for (unsigned n = 0; n < num_iters; ++n) { // find min RMS height change across N samples
			unsigned const x1(border + (rgen.rand()%xend)), y1(border + (rgen.rand()%yend));
			unsigned const x2(x1 + width), y2(y1 + height);
			if (overlaps_used (x1-slope_width, y1-slope_width, x2+slope_width, y2+slope_width)) continue; // skip if plot expanded by slope_width overlaps an existing city
			if (any_underwater(x1, y1, x2, y2, CHECK_HEIGHT_BORDER_ONLY))   continue; // skip
			float const diff(get_rms_height_diff(x1, y1, x2, y2));
			if (num_cands == 0 || diff < best_diff) {x_llc = x1; y_llc = y1; best_diff = diff;}
			if (++num_cands == num_samples) break; // done
		} // for n
		if (num_cands == 0) return 0;
		cout << "City cands: " << num_cands << ", diff: " << best_diff << ", loc: " << x_llc << "," << y_llc << endl;
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


enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, NUM_RD_TYPES};

colorRGBA const road_color(WHITE); // all road parts are the same color, to make the textures match

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

	template<typename T> static void add_flat_road_quad(T const &r, quad_batch_draw &qbd, float ar) { // z1 == z2
		float const z(r.z1()); // z1
		point const pts[4] = {point(r.x1(), r.y1(), z), point(r.x2(), r.y1(), z), point(r.x2(), r.y2(), z), point(r.x1(), r.y2(), z)};
		qbd.add_quad_pts(pts, road_color, plus_z, r.get_tex_range(ar));
	}

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
		float get_length() const {return (d[dim][1] - d[dim][0]);}
		float get_height() const {return (d[2  ][1] - d[2  ][0]);}
	};
	struct road_seg_t : public road_t {
		road_seg_t(cube_t const &c, bool dim_, bool slope_=0) : road_t(c, dim_, slope_) {}
		tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, -ar, (dim ? -1.0 : 1.0), 0, dim);}

		void add_road_quad(quad_batch_draw &qbd, float ar) const { // specialized here for sloped roads
			if (d[2][0] == d[2][1]) {add_flat_road_quad(*this, qbd, ar); return;}
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

	struct draw_state_t {
		shader_t s;
		vector3d xlate;
		bool use_smap, use_bmap;
	private:
		quad_batch_draw qbd_batched[NUM_RD_TYPES];
		bool emit_now;
		float ar;

	public:
		draw_state_t() : xlate(zero_vector), use_smap(0), use_bmap(0), emit_now(0), ar(1.0) {}
		void begin_tile(point const &pos) {emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));}

		void pre_draw() {
			if (use_smap) {
				setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 1, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
				s.add_uniform_float("z_bias", cobj_z_bias);
				s.add_uniform_float("pcf_offset", 10.0*shadow_map_pcf_offset);
			}
			ar = city_params.get_road_ar();
		}
		void post_draw() {
			emit_now = 0;
			if (use_smap) {s.end_shader();}
			setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 0, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
			
			for (unsigned i = 0; i < NUM_RD_TYPES; ++i) { // only unshadowed blocks
				road_mat_mrg.set_texture(i);
				qbd_batched[i].draw_and_clear();
			}
			s.end_shader();
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
			}
			else {qbd_batched[type_ix].add_quads(cache);} // add non-shadow blocks for drawing later
		}
	};

	class road_network_t {
		vector<road_t> roads; // full overlapping roads, for collisions, etc.
		vector<road_seg_t> segs; // non-overlapping road segments, for drawing with textures
		vector<road_isec_t> isecs[3]; // for drawing with textures: {4-way, 3-way, 2-way}
		vector<road_plot_t> plots; // plots of land that can hold buildings
		cube_t bcube;
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
		road_network_t() : bcube(all_zeros) {}
		road_network_t(cube_t const &bcube_) : bcube(bcube_) {bcube.d[2][1] += ROAD_HEIGHT;} // make it nonzero size
		cube_t const &get_bcube() const {return bcube;}
		void set_bcube(cube_t const &bcube_) {bcube = bcube_;}
		unsigned num_roads() const {return roads.size();}
		bool empty() const {return roads.empty();}

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
			float const half_width(0.5*road_width), zval(region.d[2][0] + ROAD_HEIGHT);
			float const rx1(region.x1() + half_width), rx2(region.x2() - half_width), ry1(region.y1() + half_width), ry2(region.y2() - half_width); // shrink to include centerlines
			float road_pitch(road_width + road_spacing);
			int const num_x_roads((rx2 - rx1)/road_pitch);
			road_pitch = 0.9999*(rx2 - rx1)/num_x_roads; // auto-calculate, round down slightly to avoid FP error

			// create a grid, for now; crossing roads will overlap
			for (float x = rx1; x < rx2; x += road_pitch) {
				roads.emplace_back(point(x, region.y1(), zval), point(x, region.y2(), zval), road_width, false);
			}
			unsigned const num_x(roads.size());

			for (float y = ry1; y < ry2; y += road_pitch) {
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
			p1.z   = bcube1.d[2][1];
			p2.z   = bcube2.d[2][1];
			p1[!dim] = p2[!dim] = conn_pos;
			p1[ dim] = bcube1.d[dim][ dir];
			p2[ dim] = bcube2.d[dim][!dir];
			bool const slope((p1.z < p2.z) ^ dir);
			road_t const road(p1, p2, road_width, dim, slope);
			unsigned const x1(hq.get_x_pos(road.x1())), y1(hq.get_y_pos(road.y1())), x2(hq.get_x_pos(road.x2())), y2(hq.get_y_pos(road.y2()));

			if (check_only) { // only need to do these checks in this case
				if (rn1 && !rn1->check_valid_conn_intersection(road, dim,  dir)) return -1.0; // invalid, don't make any changes
				if (rn2 && !rn2->check_valid_conn_intersection(road, dim, !dir)) return -1.0;

				for (auto b = blockers.begin(); b != blockers.end(); ++b) {
					if ((rn1 && b->contains_cube(bcube1)) || (rn2 && b->contains_cube(bcube2))) continue; // skip current cities
					// create an intersection if blocker is a road, and happens to be the same elevation?
					if (b->intersects_xy(road)) return -1.0; // bad intersection, fail
				}
				if (hq.any_underwater(x1, y1, x2+1, y2+1)) return -1.0; // underwater (Note: bounds check is done here)
			}
			// FIXME: make road follow terrain countour (could do this in split_connector_roads())?
			float const tot_dz(hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope], road.d[2][!slope], dim, check_only));
			if (check_only) return tot_dz; // return without creating the connection
			if (rn1) {rn1->insert_conn_intersection(road, dim,  dir);}
			if (rn2) {rn2->insert_conn_intersection(road, dim, !dir);}
			float const blocker_padding(max(city_params.road_spacing, 2.0f*city_params.road_border*max(DX_VAL, DY_VAL))); // use road_spacing?
			roads.push_back(road);
			blockers.push_back(road);
			blockers.back().expand_by(blocker_padding); // add extra padding
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
				float const len(r->get_length()), z1(r->d[2][slope]), z2(r->d[2][!slope]);
				assert(len > 0.0);
				unsigned const num_segs(ceil(len/road_spacing));
				cube_t c(*r); // start by copying the road's bcube
				
				for (unsigned n = 0; n < num_segs; ++n) {
					c.d[d][1] = min(r->d[d][1], (c.d[d][0] + road_spacing)); // clamp to original road end
					for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
					if (c.d[2][1] < c.d[2][0]) {swap(c.d[2][0], c.d[2][1]);} // swap zvals if needed
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
		void draw(draw_state_t &dstate) {
			if (empty()) return;
			cube_t const bcube_x(bcube + dstate.xlate);
			if (!camera_pdu.cube_visible(bcube_x)) return; // VFC
			if (!dist_less_than(camera_pdu.pos, bcube_x.closest_pt(camera_pdu.pos), get_draw_tile_dist())) return; // too far

			for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
				if (!camera_pdu.cube_visible(b->bcube + dstate.xlate)) continue; // VFC
				dstate.begin_tile(b->bcube.get_cube_center());
				dstate.draw_road_region(segs,  b->ranges[TYPE_RSEG], b->quads[TYPE_RSEG], TYPE_RSEG); // road segments
				dstate.draw_road_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PLOT], TYPE_PLOT); // plots
				for (unsigned i = 0; i < 3; ++i) {dstate.draw_road_region(isecs[i], b->ranges[TYPE_ISEC2 + i], b->quads[TYPE_ISEC2 + i], (TYPE_ISEC2 + i));} // intersections
			} // for b
		}
	}; // road_network_t

	vector<road_network_t> road_networks; // one per city
	road_network_t global_rn; // connects cities together; no plots
	draw_state_t dstate;
	rand_gen_t rgen;

	static float rgen_uniform(float val1, float val2, rand_gen_t &rgen) {return (val1 + (val2 - val1)*rgen.rand_float());}

public:
	void gen_roads(cube_t const &region, float road_width, float road_spacing) {
		timer_t timer("Gen Roads");
		road_networks.push_back(road_network_t(region));
		if (!road_networks.back().gen_road_grid(road_width, road_spacing)) {road_networks.pop_back();}
		else {cout << "Roads: " << road_networks.back().num_roads() << endl;}
	}
	bool connect_two_cities(unsigned city1, unsigned city2, vector<cube_t> &blockers, heightmap_query_t &hq, float road_width) {
		assert(city1 < road_networks.size() && city2 < road_networks.size());
		assert(city1 != city2); // check for self reference
		cout << "Connect city " << city1 << " and " << city2 << endl;
		road_network_t &rn1(road_networks[city1]), &rn2(road_networks[city2]);
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		assert(!bcube1.intersects_xy(bcube2));
		float const min_edge_dist(1.1*road_width), min_jog(2.0*road_width), half_width(0.5*road_width);
		// Note: cost function should include road length, number of jogs, total elevation change, and max slope

		for (unsigned d = 0; d < 2; ++d) { // try for single segment
			float const shared_min(max(bcube1.d[d][0], bcube2.d[d][0])), shared_max(min(bcube1.d[d][1], bcube2.d[d][1]));
			
			if (shared_max - shared_min > min_edge_dist) { // can connect with single road segment in dim !d, if the terrain in between is passable
				cout << "Shared dim " << d << endl;
				float const val1(shared_min+0.5*road_width), val2(shared_max-0.5*road_width);
				float best_conn_pos(0.0), best_cost(-1.0);

				for (unsigned n = 0; n < 20; ++n) { // make up to 20 attempts at connecting the cities with a straight line
					float const conn_pos(rgen_uniform(val1, val2, rgen)); // chose a random new connection point and try it
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, hq, road_width, conn_pos, !d, 1)); // check_only=1
					if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost;}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					cout << "Single segment cost: " << best_cost << endl;
					global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, hq, road_width, best_conn_pos, !d, 0); // check_only=0; actually make the change
					return 1;
				}
			}
		} // for d
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		bool const dx(center1.x < center2.x), dy(center1.y < center2.y), slope(dx ^ dy); // dx=1, dy=1
		cube_t const bc[2] = {bcube1, bcube2};
		
		// FIXME: make this work with partial overlap case
		if ((bc[dx].x1() - bc[!dx].x2() > min_jog) && (bc[dy].y1() - bc[!dy].y2() > min_jog)) { // connect with two road segments using a jog
			bool const inv_dim(rgen.rand()&1);
			cout << "Try connect using jog in dim " << inv_dim << endl;

			for (unsigned d = 0; d < 2; ++d) { // x-then-y vs. y-then-x
				bool const fdim((d != 0) ^ inv_dim); // first segment dim: 0=x first, 1=y first (== 0)
				bool const range_dir1(fdim ? dy : dx), range_dir2(fdim ? dx : dy); // 1/1
				cube_t region1(bcube1), region2(bcube2);
				region1.d[ fdim][!range_dir1] = region1.d[ fdim][ range_dir1];
				region2.d[!fdim][ range_dir2] = region2.d[!fdim][!range_dir2];
				region1.d[!fdim][0] += min_edge_dist; region1.d[!fdim][1] -= min_edge_dist;
				region2.d[ fdim][0] += min_edge_dist; region1.d[ fdim][1] -= min_edge_dist;
				float const xmin(region1.d[!fdim][0]), xmax(region1.d[!fdim][1]), ymin(region2.d[fdim][0]), ymax(region2.d[fdim][1]);
				float best_xval(0.0), best_yval(0.0), best_cost(-1.0);
				cube_t best_int_cube;

				for (unsigned n = 0; n < 20; ++n) { // make up to 20 attempts at connecting the cities with a single jog
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
					cout << "Double segment cost: " << best_cost << endl;
					//cout << TXT(best_xval) << TXT(best_yval) << TXT(fdim) << ", int_cube: " << best_int_cube.str() << endl;
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), hq); // do this first to improve flattening
					global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, hq, road_width, (fdim ? best_xval : best_yval),  fdim, 0); // check_only=0
					global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2, hq, road_width, (fdim ? best_yval : best_xval), !fdim, 0); // check_only=0
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
				cout << "Trying to connect city " << i << " to city " << j << ": " << success << endl;
			} // for j
		} // for i
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
		dstate.use_smap = shadow_map_enabled();
		dstate.xlate    = xlate;
		dstate.pre_draw();
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate);}
		global_rn.draw(dstate);
		dstate.post_draw();
		glDepthFunc(GL_LESS);
		fgPopMatrix();
	}
}; // city_road_gen_t


class city_gen_t : public city_plot_gen_t {

	city_road_gen_t road_gen;

public:
	bool gen_city(city_params_t const &params, cube_t &cities_bcube) {
		timer_t t("Choose City Location");
		unsigned x1(0), y1(0);
		if (!find_best_city_location(params.city_size, params.city_size, params.city_border, params.slope_width, params.num_samples, x1, y1)) return 0;
		unsigned const x2(x1 + params.city_size), y2(y1 + params.city_size);
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		if (params.roads_enabled()) {road_gen.gen_roads(pos_range, params.road_width, params.road_spacing);}
		return 1;
	}
	void gen_cities(city_params_t const &params) {
		if (params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params, cities_bcube);}
		bool const is_const_zval(cities_bcube.d[2][0] == cities_bcube.d[2][1]);
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube, is_const_zval);}
		road_gen.connect_all_cities(heightmap, xsize, ysize, params.road_width, params.road_spacing);
		road_gen.gen_tile_blocks();
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_road_bcubes(bcubes);}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_plot_bcubes(bcubes);}

	bool check_city_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (check_plot_sphere_coll(pos, radius, xy_only)) return 1;
		return road_gen.check_road_sphere_coll(pos, radius, xy_only);
	}
	void draw(bool shadow_only, int reflection_pass, vector3d const &xlate) { // for now, there are only roads
		if (!shadow_only && reflection_pass == 0) {road_gen.draw(xlate);} // roads don't cast shadows and aren't reflected in water
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

