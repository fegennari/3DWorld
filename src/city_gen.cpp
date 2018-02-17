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
float const ROAD_HEIGHT = 0.001;


extern int rand_gen_index, display_mode;
extern float water_plane_z, shadow_map_pcf_offset, cobj_z_bias;


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
	vector<cube_t> plots; // same size as used

	bool is_valid_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		return (x1 < x2 && y1 < y2 && x2 <= xsize && y2 <= ysize);
	}
	bool overlaps_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		rect_t const cur(x1, y1, x2, y2);
		for (vector<rect_t>::const_iterator i = used.begin(); i != used.end(); ++i) {if (i->has_overlap(cur)) return 1;} // simple linear iteration
		return 0;
	}
	cube_t add_plot(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) {
		cube_t bcube;
		int const dx(-int(xsize)/2), dy(-int(ysize)/2); // convert from center to LLC
		bcube.d[0][0] = get_xval(x1 + dx);
		bcube.d[0][1] = get_xval(x2 + dx);
		bcube.d[1][0] = get_yval(y1 + dy);
		bcube.d[1][1] = get_yval(y2 + dy);
		bcube.d[2][0] = bcube.d[2][1] = elevation;
		plots.push_back(bcube);
		used.emplace_back(x1, y1, x2, y2);
		return bcube;
	}
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
	vector3d const get_query_xlate() const {
		return vector3d((world_mode == WMODE_INF_TERRAIN) ? vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0) : zero_vector);
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
		return 1; // success
	}
	float flatten_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float const *const height=nullptr) {
		assert(is_valid_region(x1, y1, x2, y2));
		float const delta_h = 0.0; // for debugging in map view
		float const elevation(height ? *height : (get_avg_height(x1, y1, x2, y2) + delta_h));

		for (unsigned y = max((int)y1-(int)slope_width, 0); y < min(y2+slope_width, ysize); ++y) {
			for (unsigned x = max((int)x1-(int)slope_width, 0); x < min(x2+slope_width, xsize); ++x) {
				float const dx(max(0, max(((int)x1 - (int)x), ((int)x - (int)x2 + 1))));
				float const dy(max(0, max(((int)y1 - (int)y), ((int)y - (int)y2 + 1))));
				float const mix(sqrt(dx*dx + dy*dy)/slope_width);
				float &h(heightmap[y*xsize + x]);
				h = mix*h + (1.0 - mix)*elevation;
			}
		}
		return elevation;
	}
	bool check_plot_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (plots.empty()) return 0;
		point const sc(pos - get_query_xlate());

		for (auto i = plots.begin(); i != plots.end(); ++i) {
			if (xy_only ? sphere_cube_intersect_xy(sc, radius, *i) : sphere_cube_intersect(sc, radius, *i)) return 1;
		}
		return 0;
	}
}; // city_plot_gen_t


struct road_mat_t {
	unsigned tid;
	colorRGBA color;
	road_mat_t(unsigned tid_=0, colorRGBA const &color_=WHITE) : tid(tid_), color(color_) {}
};

road_mat_t const road_mats[4] = {
	road_mat_t(WHITE_TEX, colorRGBA(0.06, 0.06, 0.06, 1.0)), // xsegs
	road_mat_t(WHITE_TEX, colorRGBA(0.06, 0.06, 0.06, 1.0)), // ysegs
	road_mat_t(TILE_TEX, colorRGBA(0.06, 0.06, 0.06, 1.0)), // isecs
	road_mat_t(STUCCO_TEX, colorRGBA(0.12, 0.12, 0.12, 1.0)), // plots
};


class city_road_gen_t {

	struct draw_state_t {
		shader_t s;
		vector3d xlate;
		bool use_smap, use_bmap;
	private:
		quad_batch_draw qbd_batched[4], qbd_shadow[4];
		unsigned type_ix;
		bool emit_now;

		quad_batch_draw &get_qbd() {
			assert(type_ix < 4);
			return (emit_now ? qbd_shadow[type_ix] : qbd_batched[type_ix]);
		}
		void draw_cur_block() {
			assert(type_ix < 4);
			select_texture(road_mats[type_ix].tid);
			get_qbd().draw_and_clear();
		}
	public:
		draw_state_t() : xlate(zero_vector), use_smap(0), use_bmap(0), type_ix(0), emit_now(0) {}
		void begin_tile(point const &pos) {emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));}
		void begin_block(unsigned type_ix_) {type_ix = type_ix_;}

		void pre_draw() {
			if (use_smap) {
				setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 1, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
				s.add_uniform_float("z_bias", cobj_z_bias);
				s.add_uniform_float("pcf_offset", 10.0*shadow_map_pcf_offset);
			}
		}
		void add_cube(cube_t const &c, tex_range_t const &tr) {
			float const z(0.5*(c.d[2][0] + c.d[2][1]));
			point const pts[4] = {point(c.d[0][0], c.d[1][0], z), point(c.d[0][1], c.d[1][0], z), point(c.d[0][1], c.d[1][1], z), point(c.d[0][0], c.d[1][1], z)};
			get_qbd().add_quad_pts(pts, road_mats[type_ix].color, plus_z, tr);
		}
		void end_block() {
			if (emit_now) {draw_cur_block();} // only shadow blocks
		}
		void post_draw() {
			emit_now = 0;
			if (use_smap) {s.end_shader();}
			setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 0, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
			for (unsigned i = 0; i < 4; ++i) {begin_block(i); draw_cur_block();} // only unshadowed blocks
			s.end_shader();
		}
	};

	struct road_t {
		cube_t bcube;

		road_t() {}
		road_t(float x1, float x2, float y1, float y2, float z1, float z2) : bcube(x1, x2, y1, y2, z1, z2) {}
		road_t(point const &s, point const &e, float width) {
			assert(s != e);
			assert(width > 0.0);
			vector3d const dw(0.5*width*cross_product((e - s), plus_z).get_norm());
			point const pts[4] = {(s - dw), (s + dw), (e + dw), (e - dw)};
			bcube = cube_t(pts, 4);
		}
	};

	class road_network_t {
		vector<road_t> roads; // full overlapping roads, for collisions, etc.
		vector<cube_t> xsegs, ysegs; // non-overlapping road segments, for drawing with textures
		vector<cube_t> isecs; // for drawing with textures
		vector<cube_t> plots; // plots of land that can hold buildings
		cube_t bcube;

		static uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}

		struct cmp_by_tile { // not the most efficient solution, but no memory overhead
			bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
		};

		struct range_pair_t {
			unsigned s, e; // Note: e is one past the end
			range_pair_t() : s(0), e(0) {}
			void update(unsigned v) {
				if (s == 0 && e == 0) {s = v;} // first insert
				else {assert(s < e && v >= e);} // v must strictly increase
				e = v+1; // one past the end
			}
		};
		struct tile_block_t { // collection of road parts for a given tile
			range_pair_t ranges[4]; // {xsegs, ysegs, isecs, plots}
			cube_t bcube;
			tile_block_t(cube_t const &bcube_) : bcube(bcube_) {}
		};
		vector<tile_block_t> tile_blocks;

		void add_tile_blocks(vector<cube_t> &v, map<uint64_t, unsigned> &tile_to_block_map, unsigned type_ix) {
			assert(type_ix < 4);
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
		void gen_tile_blocks() {
			tile_blocks.clear(); // should already be empty?
			map<uint64_t, unsigned> tile_to_block_map;
			add_tile_blocks(xsegs, tile_to_block_map, 0);
			add_tile_blocks(ysegs, tile_to_block_map, 1);
			add_tile_blocks(isecs, tile_to_block_map, 2);
			add_tile_blocks(plots, tile_to_block_map, 3);
			//cout << "tile_to_block_map: " << tile_to_block_map.size() << ", tile_blocks: " << tile_blocks.size() << endl;
		}

	public:
		road_network_t(cube_t const &bcube_) : bcube(bcube_) {bcube.d[2][1] += ROAD_HEIGHT;} // make it nonzero size
		unsigned num_roads() const {return roads.size();}
		void clear() {roads.clear(); xsegs.clear(); ysegs.clear(); isecs.clear(); plots.clear(); tile_blocks.clear();}

		bool gen_road_grid(cube_t const &region, float road_width, float road_spacing) {
			vector3d const size(region.get_size());
			assert(size.x > 0.0 && size.y > 0.0);
			float const half_width(0.5*road_width), road_pitch(road_width + road_spacing);
			float const zval(region.d[2][0] + ROAD_HEIGHT);

			// create a grid, for now; crossing roads will overlap
			for (float x = region.d[0][0]+half_width; x < region.d[0][1]-half_width; x += road_pitch) { // shrink to include centerlines
				roads.emplace_back(point(x, region.d[1][0], zval), point(x, region.d[1][1], zval), road_width);
			}
			unsigned const num_x(roads.size());

			for (float y = region.d[1][0]+half_width; y < region.d[1][1]-half_width; y += road_pitch) { // shrink to include centerlines
				roads.emplace_back(point(region.d[0][0], y, zval), point(region.d[0][1], y, zval), road_width);
			}
			unsigned const num_y(roads.size() - num_x);
			if (num_x <= 1 || num_y <= 1) {clear(); return 0;} // not enough space for roads

			// create road segments and intersections
			xsegs.reserve(num_x*(num_y-1));
			ysegs.reserve((num_x-1)*num_y);
			isecs.reserve(num_x*num_y);
			plots.reserve((num_x-1)*(num_y-1));

			for (unsigned x = 0; x < num_x; ++x) {
				for (unsigned y = num_x; y < roads.size(); ++y) {
					cube_t const &rx(roads[x].bcube), &ry(roads[y].bcube);
					isecs.emplace_back(rx.d[0][0], rx.d[0][1], ry.d[1][0], ry.d[1][1], zval, zval); // intersections
					if (x+1 < num_x) { // skip last y segment
						cube_t const &rxn(roads[x+1].bcube);
						ysegs.emplace_back(rx.d[0][1], rxn.d[0][0], ry.d[1][0], ry .d[1][1], zval, zval); // y-segments
					}
					if (y+1 < roads.size()) { // skip last x segment
						cube_t const &ryn(roads[y+1].bcube);
						xsegs.emplace_back(rx.d[0][0], rx .d[0][1], ry.d[1][1], ryn.d[1][0], zval, zval); // x-segments

						if (x+1 < num_x) { // skip last y segment
							cube_t const &rxn(roads[x+1].bcube);
							plots.emplace_back(rx.d[0][1], rxn.d[0][0], ry.d[1][1], ryn.d[1][0], zval, zval); // plots between roads
						}
					}
				} // for y
			} // for x
			gen_tile_blocks();
			return 1;
		}
		void get_road_bcubes(vector<cube_t> &bcubes) const {
			for (auto r = roads.begin(); r != roads.end(); ++r) {bcubes.push_back(r->bcube);}
		}
		void get_plot_bcubes(vector<cube_t> &bcubes) const {
			for (auto r = plots.begin(); r != plots.end(); ++r) {bcubes.push_back(*r);}
		}
		void draw(draw_state_t &dstate) const {
			cube_t const bcube_x(bcube + dstate.xlate);
			if (!camera_pdu.cube_visible(bcube_x)) return; // VFC
			if (!dist_less_than(camera_pdu.pos, bcube_x.closest_pt(camera_pdu.pos), get_draw_tile_dist())) return; // too far
			colorRGBA const road_color(0.06, 0.06, 0.06, 1.0), int_color(0.06, 0.06, 0.06, 1.0), plot_color(0.12, 0.12, 0.12, 1.0);

			for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
				if (!camera_pdu.cube_visible(b->bcube + dstate.xlate)) continue; // VFC
				dstate.begin_tile(b->bcube.get_cube_center());
				draw_road_region(xsegs, b->ranges[0], dstate, road_color, 0); // x road segments
				draw_road_region(ysegs, b->ranges[1], dstate, road_color, 1); // y road segments
				draw_road_region(isecs, b->ranges[2], dstate, int_color,  2); // intersections
				draw_road_region(plots, b->ranges[3], dstate, plot_color, 3); // plots
			} // for b
		}
		private:
			void draw_road_region(vector<cube_t> const &v, range_pair_t const &rp, draw_state_t &dstate, colorRGBA const &color, unsigned type_ix) const {
				assert(rp.s <= rp.e && rp.e <= v.size());
				dstate.begin_block(type_ix);
				for (unsigned i = rp.s; i < rp.e; ++i) {dstate.add_cube(v[i], tex_range_t());}
				dstate.end_block();
			}
	};

	vector<road_network_t> road_networks;
	draw_state_t dstate;

public:
	void gen_roads(cube_t const &region, float road_width, float road_spacing) {
		timer_t timer("Gen Roads");
		road_networks.push_back(road_network_t(region));
		if (!road_networks.back().gen_road_grid(region, road_width, road_spacing)) {road_networks.pop_back();}
		else {cout << "Roads: " << road_networks.back().num_roads() << endl;}
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_road_bcubes(bcubes);}
	}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const {
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_plot_bcubes(bcubes);}
	}
	void draw(vector3d const &xlate) { // non-const because qbd is modified
		if (road_networks.empty()) return;
		//timer_t timer("Draw Roads");
		fgPushMatrix();
		translate_to(xlate);
		glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
		dstate.use_smap = shadow_map_enabled();
		dstate.xlate    = xlate;
		select_texture(WHITE_TEX); // FIXME: ROAD_TEX
		dstate.pre_draw();
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate);}
		dstate.post_draw();
		glDepthFunc(GL_LESS);
		fgPopMatrix();
	}
}; // city_road_gen_t


struct city_params_t {

	unsigned num_cities, num_samples, city_size, city_border, slope_width;
	float road_width, road_spacing;

	city_params_t() : num_cities(0), num_samples(100), city_size(0), city_border(0), slope_width(0), road_width(0.0), road_spacing(0.0) {}
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


class city_gen_t : public city_plot_gen_t {

	city_road_gen_t road_gen;

public:
	bool gen_city(city_params_t const &params, cube_t &cities_bcube) {
		timer_t t("Choose City Location");
		unsigned x1(0), y1(0);
		if (!find_best_city_location(params.city_size, params.city_size, params.city_border, params.num_samples, x1, y1)) return 0;
		unsigned const x2(x1 + params.city_size), y2(y1 + params.city_size);
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		if (params.road_width > 0.0 && params.road_spacing > 0.0) {road_gen.gen_roads(pos_range, params.road_width, params.road_spacing);}
		return 1;
	}
	void gen_cities(city_params_t const &params) {
		if (params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params, cities_bcube);}
		bool const is_const_zval(cities_bcube.d[2][0] == cities_bcube.d[2][1]);
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube, is_const_zval);}
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_road_bcubes(bcubes);}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_plot_bcubes(bcubes);}

	void draw(bool shadow_only, int reflection_pass, vector3d const &xlate) { // for now, there are only roads
		if (!shadow_only && reflection_pass == 0) {road_gen.draw(xlate);} // roads don't cast shadows and aren't reflected in water
		// buildings are drawn through draw_buildings()
	}
}; // city_gen_t


city_params_t city_params;
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
	return city_gen.check_plot_sphere_coll(center, radius);
}
bool check_valid_scenery_pos(point const &pos, float radius) {
	if (check_buildings_sphere_coll(pos, radius, 1, 1)) return 0; // apply_tt_xlate=1, xy_only=1
	if (check_city_sphere_coll(pos, radius)) return 0;
	return 1;
}

