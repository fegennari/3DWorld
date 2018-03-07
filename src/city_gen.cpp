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
float const CAR_LANE_OFFSET         = 0.15; // in units of road width
vector3d const CAR_SIZE(0.28, 0.13, 0.07); // {length, width, height} in units of road width
colorRGBA const road_color          = WHITE; // all road parts are the same color, to make the textures match
unsigned  const CONN_CITY_IX((1<<16)-1); // uint16_t max

enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, NUM_RD_TYPES};
enum {TURN_NONE=0, TURN_LEFT, TURN_RIGHT};
unsigned const CONN_TYPE_NONE = 0;


extern int rand_gen_index, display_mode, animate2;
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

bool is_isect(unsigned type) {return (type >= TYPE_ISEC2 && type <= TYPE_ISEC4);}
int encode_neg_ix(unsigned ix) {return -(int(ix)+1);}
unsigned decode_neg_ix(int ix) {assert(ix < 0); return -(ix+1);}

class city_road_gen_t;

struct car_t {
	cube_t bcube, prev_bcube;
	bool dim, dir, stopped_at_light;
	unsigned char cur_road_type, color_id, turn_dir;
	unsigned short cur_city, cur_road, cur_seg;
	float dz, height, cur_speed, max_speed;

	car_t() : bcube(all_zeros), dim(0), dir(0), stopped_at_light(0), cur_road_type(0), color_id(0), turn_dir(TURN_NONE),
		cur_city(0), cur_road(0), cur_seg(0), dz(0.0), height(0.0), cur_speed(0.0), max_speed(0.0) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	point get_center() const {return bcube.get_cube_center();}
	unsigned get_orient() const {return (2*dim + dir);}
	float get_length() const {return (bcube.d[dim][1] - bcube.d[dim][0]);}
	bool is_stopped () const {return (cur_speed == 0.0);}
	bool is_parked  () const {return (max_speed == 0.0);}
	void park() {cur_speed = max_speed = 0.0;}

	bool operator<(car_t const &c) const { // sort spatially for collision detection and drawing
		if (cur_city != c.cur_city) return (cur_city < c.cur_city);
		if (cur_road != c.cur_road) return (cur_road < c.cur_road);
		if (cur_road_type != c.cur_road_type) return (cur_road_type < c.cur_road_type);
		if (cur_seg  != c.cur_seg ) return (cur_seg  < c.cur_seg );
		return (bcube.get_cube_center() < c.bcube.get_cube_center());
	}
	string str() const {
		std::ostringstream oss;
		oss << "Car " << TXT(dim) << TXT(dir) << TXT(cur_city) << TXT(cur_road) << TXT(cur_seg) << TXT(dz) << TXT(max_speed) << TXT(cur_speed)
			<< "cur_road_type=" << unsigned(cur_road_type) << " color=" << unsigned(color_id) << " bcube=" << bcube.str();
		return oss.str();
	}
	void move(float speed_mult) {
		prev_bcube = bcube;
		if (cur_speed == 0.0) return;
		assert(speed_mult >= 0.0 && cur_speed > 0.0 && cur_speed <= max_speed);
		float dist(cur_speed*speed_mult);
		min_eq(dist, 0.25f*city_params.road_width); // limit to half a car length to prevent cars from crossing an intersection in a single frame
		move_by((dir ? 1.0 : -1.0)*dist);
	}
	void accelerate(float mult=0.02) {cur_speed = min(max_speed, (cur_speed + mult*fticks*max_speed));}
	void decelerate(float mult=0.05) {cur_speed = max(0.0f, (cur_speed - mult*fticks*max_speed));}
	void decelerate_fast() {decelerate(10.0);} // Note: large decel to avoid stopping in an intersection
	void move_by(float val) {bcube.d[dim][0] += val; bcube.d[dim][1] += val;}
	bool check_collision(car_t &c, city_road_gen_t const &road_gen);
};

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
	unsigned road_ix;
	bool dim; // dim the road runs in
	bool slope; // 0: z1 applies to first (lower) point; 1: z1 applies to second (upper) point

	road_t(cube_t const &c, bool dim_, bool slope_=0, unsigned road_ix_=0) : cube_t(c), road_ix(road_ix_), dim(dim_), slope(slope_) {}
	road_t(point const &s, point const &e, float width, bool dim_, bool slope_=0, unsigned road_ix_=0) : road_ix(road_ix_), dim(dim_), slope(slope_) {
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
	unsigned short road_ix, conn_ix[2], conn_type[2]; // {dim=0, dim=1}
	void init_ixs() {conn_ix[0] = conn_ix[1] = 0; conn_type[0] = conn_type[1] = CONN_TYPE_NONE;}
	road_seg_t(road_t const &r, unsigned rix) : road_t(r), road_ix(rix) {init_ixs();}
	road_seg_t(cube_t const &c, unsigned rix, bool dim_, bool slope_=0) : road_t(c, dim_, slope_), road_ix(rix) {init_ixs();}
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

namespace stoplight_ns {

	enum {RED=0, GREEN=1, YELLOW=2}; // colors, unused (only have stop and go states anyway)
	enum {EGL=0, EGWG, WGL, NGL, NGSG, SGL, NUM_STATE}; // E=car moving east, W=west, N=sorth, S=south, G=straight|right, L=left turn
	float const state_times[NUM_STATE] = {4.0, 6.0, 4.0, 4.0, 6.0, 4.0}; // in seconds
	unsigned const st_r_orient_masks[NUM_STATE] = {2, 3, 1, 8, 12, 4}; // {W=1, E=2, S=4, N=8}, for straight and right turns
	unsigned const left_orient_masks[NUM_STATE] = {2, 0, 1, 8, 0,  4}; // {W=1, E=2, S=4, N=8}, for left turns only

	rand_gen_t stoplight_rgen;

	class stoplight_t {
		uint8_t num_conn, conn, cur_state;
		float cur_state_ticks;

		void next_state() {
			++cur_state;
			if (cur_state == NUM_STATE) {cur_state = 0;} // wraparound
		}
		void advance_state() {
			if (num_conn == 4) {next_state();} // all states are valid for 4-way intersections
			else { // 3-way intersection
				assert(num_conn == 3);
				while (1) {
					next_state();
					bool valid(0);
					switch (conn) {
					case 7 : {bool const allow[6] = {0,1,1,1,0,0}; valid = allow[cur_state]; break;} // no +y
					case 11: {bool const allow[6] = {1,1,0,0,0,1}; valid = allow[cur_state]; break;} // no -y
					case 13: {bool const allow[6] = {1,0,0,1,1,0}; valid = allow[cur_state]; break;} // no +x
					case 14: {bool const allow[6] = {0,0,1,0,1,1}; valid = allow[cur_state]; break;} // no -x
					default: assert(0);
					}
					if (valid) break;
				} // end while
			}
			cur_state_ticks = 0.0; // reset for this state
		}
		void run_update_logic() {
			assert(cur_state < NUM_STATE);
			if (cur_state_ticks < TICKS_PER_SECOND*state_times[cur_state]) return; // keep existing state
			advance_state();
		}
	public:
		stoplight_t() : num_conn(0), conn(0), cur_state(0), cur_state_ticks(0.0) {}
		void init(uint8_t num_conn_, uint8_t conn_) {
			num_conn = num_conn_; conn = conn_;
			if (num_conn == 2) return; // nothing else to do
			cur_state = stoplight_rgen.rand() % NUM_STATE; // start at a random state
			advance_state(); // make sure cur_state is valid
			cur_state_ticks = TICKS_PER_SECOND*state_times[cur_state]*stoplight_rgen.rand_float(); // start at a random time within this state
		}
		void next_frame() {
			if (num_conn == 2) return; // nothing else to do
			cur_state_ticks += fticks;
			run_update_logic();
		}
		bool red_light(bool dim, bool dir, unsigned turn) const {
			assert(cur_state < NUM_STATE);
			assert(turn == TURN_LEFT || turn == TURN_RIGHT || turn == TURN_NONE);
			if (num_conn == 2) return 0; // 2-way intersection, no cross traffic
			unsigned const mask(((turn == TURN_LEFT) ? left_orient_masks : st_r_orient_masks)[cur_state]);
			unsigned const orient(2*dim + dir); // {W, E, S, N}
			return (((1<<orient) & mask) == 0);
		}
		bool red_or_yellow_light(bool dim, bool dir, unsigned turn) const {
			if (red_light(dim, dir, turn)) return 1;
			if (num_conn == 2) return 0;
			stoplight_t future_self(*this);
			float const yellow_light_time(2.0);
			future_self.cur_state_ticks += TICKS_PER_SECOND*yellow_light_time;
			future_self.run_update_logic();
			return future_self.red_light(dim, dir, turn);
		}
	};
} // end stoplight_ns

struct road_isec_t : public cube_t {
	uint8_t num_conn, conn; // connected roads in {-x, +x, -y, +y}
	short rix_xy[2], conn_ix[4]; // pos=cur city road, neg=global road; always segment ix
	stoplight_ns::stoplight_t stoplight; // Note: not always needed, maybe should be by pointer/index?

	road_isec_t(cube_t const &c, int rx, int ry, uint8_t conn_) : cube_t(c), conn(conn_) {
		rix_xy[0] = rx; rix_xy[1] = ry; conn_ix[0] = conn_ix[1] = conn_ix[2] = conn_ix[3] = 0;
		if (conn == 15) {num_conn = 4;} // 4-way
		else if (conn == 7 || conn == 11 || conn == 13 || conn == 14) {num_conn = 3;} // 3-way
		else if (conn == 5 || conn == 6  || conn == 9  || conn == 10) {num_conn = 2;} // 2-way
		else {assert(0);}
		stoplight.init(num_conn, conn);
	}
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
	void make_4way() {num_conn = 4; conn = 15;}
	void next_frame() {stoplight.next_frame();}
	bool red_light(car_t const &car) const {return stoplight.red_light(car.dim, car.dir, car.turn_dir);}
	bool red_or_yellow_light(car_t const &car) const {return stoplight.red_or_yellow_light(car.dim, car.dir, car.turn_dir);}
};

struct road_plot_t : public cube_t {
	road_plot_t(cube_t const &c) : cube_t(c) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, ar, ar);}
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
		unsigned city_id, cluster_id;
		//string city_name; // future work

		// use only for the global road network
		struct city_id_pair_t {
			unsigned id[2]; // lo, hi
			city_id_pair_t(unsigned c1, unsigned c2) {id[0] = c1; id[1] = c2;}
		};
		vector<city_id_pair_t> road_to_city; // indexed by road ID
		vector<vector<unsigned>> city_to_seg; // maps city_id to set of road segments connecting to that city

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
		road_network_t() : bcube(all_zeros), city_id(CONN_CITY_IX), cluster_id(0) {} // global road network ctor
		road_network_t(cube_t const &bcube_, unsigned city_id_) : bcube(bcube_), city_id(city_id_), cluster_id(0) {bcube.d[2][1] += ROAD_HEIGHT;} // make it nonzero size
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
				roads.emplace_back(point(x, region.y1(), zval), point(x, region.y2(), zval), road_width, true);
			}
			unsigned const num_x(roads.size());

			for (float y = ry1; y < ry2; y += road_pitch_y) {
				roads.emplace_back(point(region.x1(), y, zval), point(region.x2(), y, zval), road_width, false);
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
					isecs[num_conn - 2].emplace_back(cube_t(rx.x1(), rx.x2(), ry.y1(), ry.y2(), zval, zval), y, x, conn); // intersections
					
					if (!LX) { // skip last y segment
						cube_t const &rxn(roads[x+1]);
						segs.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y1(), ry.y2(), zval, zval), y, false); // y-segments
					}
					if (!LY) { // skip last x segment
						cube_t const &ryn(roads[y+1]);
						segs.emplace_back(cube_t(rx.x1(), rx.x2(), ry.y2(), ryn.y1(), zval, zval), x, true); // x-segments

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
			isec.make_4way(); // all connected
			isecs[2].push_back(isec); // add as 4-way intersection
			isecs[1][int3_ix] = isecs[1].back(); // remove original 3-way intersection
			isecs[1].pop_back();
		}

		struct road_ixs_t {
			vector<unsigned> seg_ixs, isec_ixs[3][2]; // {2-way, 3-way, 4-way} x {X, Y}
		};
		template<typename T> int search_for_adj(vector<T> const &v, vector<unsigned> const &ixs, cube_t const &bcube, bool dim, bool dir, bool debug=0) const {
			if (debug) {cout << "bcube: " << bcube.str() << " " << TXT(dim) << TXT(dir) << TXT(ixs.size()) << endl;}
			for (auto i = ixs.begin(); i != ixs.end(); ++i) {
				assert(*i < v.size());
				cube_t const &c(v[*i]);
				if (debug) {cout << "c[" << *i << "]: " << c.str() << " " << (c.d[dim][!dir]==bcube.d[dim][dir]) << (c.d[!dim][0]==bcube.d[!dim][0]) << (c.d[!dim][1]==bcube.d[!dim][1]) << endl;}
				if (c.d[dim][!dir] != bcube.d[dim][dir]) continue; // no shared edge
				if (c.d[!dim][0] != bcube.d[!dim][0] || c.d[!dim][1] != bcube.d[!dim][1]) continue; // no shared edge in other dim
				return *i; // there can be only one
			} // for i
			return -1; // not found
		}
		void calc_ix_values(road_network_t const &global_rn) {
			// now that the segments and intersections are in order, we can fill in the IDs; first, create a mapping from road to segments and intersections
			bool const is_global_rn(&global_rn == this);
			assert(road_to_city.size() == (is_global_rn ? roads.size() : 0));
			vector<road_ixs_t> by_ix(roads.size()); // maps road_ix to list of seg_ix values
			
			if (is_global_rn) {
				unsigned num_cities(0);
				for (auto r = road_to_city.begin(); r != road_to_city.end(); ++r) {
					for (unsigned d = 0; d < 2; ++d) {if (r->id[d] != CONN_CITY_IX) {max_eq(num_cities, r->id[d]+1);}}
				}
				city_to_seg.resize(num_cities);
			}
			for (unsigned i = 0; i < segs.size(); ++i) {
				unsigned const ix(segs[i].road_ix);
				assert(ix < by_ix.size());
				assert(roads[ix].dim == segs[i].dim);
				by_ix[ix].seg_ixs.push_back(i);
			}
			for (unsigned n = 0; n < 3; ++n) { // {2-way, 3-way, 4-way}
				for (unsigned i = 0; i < isecs[n].size(); ++i) {
					for (unsigned d = 0; d < 2; ++d) { // {x, y}
						int ix(isecs[n][i].rix_xy[d]);
						//cout << TXT(n) << TXT(i) << TXT(d) << TXT(ix) << endl;
						
						if (ix < 0) { // global connector road
							ix = decode_neg_ix(ix);
							assert((unsigned)ix < global_rn.roads.size());
							continue; // connector road, nothing else to do here?
						}
						else {
							assert((unsigned)ix < roads.size());
							assert(roads[ix].dim == (d != 0));
							by_ix[ix].isec_ixs[n][d].push_back(i);
						}
					} // for d
				} // for i
			} // for n

			// next, connect segments and intersections together by index, using roads as an acceleration structure
			for (unsigned i = 0; i < segs.size(); ++i) {
				road_seg_t &seg(segs[i]);
				road_ixs_t const &rix(by_ix[seg.road_ix]);

				for (unsigned dir = 0; dir < 2; ++dir) { // dir
					bool found(0);
					int const seg_ix(search_for_adj(segs, rix.seg_ixs, seg, seg.dim, (dir != 0)));
					//cout << TXT(i) << TXT(seg.road_ix) << TXT(dir) << TXT(seg.dim) << TXT(seg_ix) << "sz: " << rix.seg_ixs.size() << endl;
					if (seg_ix >= 0) {assert(seg_ix != i); seg.conn_ix[dir] = seg_ix; seg.conn_type[dir] = TYPE_RSEG; found = 1;} // found segment

					for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
						int const isec_ix(search_for_adj(isecs[n], rix.isec_ixs[n][seg.dim], seg, seg.dim, (dir != 0)));
						//cout << TXT(n) << TXT(isec_ix) << "sz: " << rix.isec_ixs[n][seg.dim].size() << endl;
						if (isec_ix >= 0) {assert(!found); seg.conn_ix[dir] = isec_ix; seg.conn_type[dir] = (TYPE_ISEC2 + n); found = 1;} // found intersection
					}
					if (is_global_rn && !found) { // connection to a city
						seg.conn_type[dir] = TYPE_ISEC3; // always connects to a 3-way intersection within the city
						assert(seg.road_ix < road_to_city.size());
						unsigned const city(road_to_city[seg.road_ix].id[dir]);
						assert(city != CONN_CITY_IX); // internal segments should be connected and not get here
						seg.conn_ix[dir] = city; // hacked in as city_id, need to handle on the other end (or is it unused?)
						//cout << TXT(city) << TXT(i) << endl;
						assert(city < city_to_seg.size());
						city_to_seg[city].push_back(i); // add segment ID
						continue;
					}
					assert(found);
				} // for dir
			} // for i
			for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
				for (unsigned i = 0; i < isecs[n].size(); ++i) {
					road_isec_t &isec(isecs[n][i]);

					for (unsigned d = 0; d < 4; ++d) { // {-x, +x, -y, +y}
						if (!(isec.conn & (1<<d))) continue; // no connection in this position
						unsigned const dim(d>>1), dir(d&1);
						int const ix(isec.rix_xy[dim]);
						
						if (ix < 0) { // global connector road
							//cout << TXT(ix) << TXT(city_id) << endl;
							vector<unsigned> const &seg_ids(global_rn.get_segs_connecting_to_city(city_id));
							assert(!seg_ids.empty());
							int const seg_ix(global_rn.search_for_adj(global_rn.segs, seg_ids, isec, (dim != 0), (dir != 0))); // global conn segment
							assert(seg_ix >= 0); // must be found
							isec.conn_ix[d] = encode_neg_ix(seg_ix);
						}
						else { // local segment
							int const seg_ix(search_for_adj(segs, by_ix[ix].seg_ixs, isec, (dim != 0), (dir != 0))); // always connects to a road segment
							//cout << TXT(n) << TXT(i) << TXT(d) << TXT(dim) << TXT(dir) << TXT(ix) << TXT(seg_ix) << "sz " << by_ix[ix].seg_ixs.size() << endl;
							assert(seg_ix >= 0); // must be found
							isec.conn_ix[d] = seg_ix;
						}
					} // for d
				} // for i
			} // for n
		}
		vector<unsigned> const &get_segs_connecting_to_city(unsigned city) const {
			assert(city < city_to_seg.size());
			return city_to_seg[city];
		}
	public:
		bool check_valid_conn_intersection(cube_t const &c, bool dim, bool dir) const {return (find_conn_int_seg(c, dim, dir) >= 0);}
		void insert_conn_intersection(cube_t const &c, bool dim, bool dir, unsigned grn_rix) { // Note: dim is the dimension of the connector road
			int const seg_id(find_conn_int_seg(c, dim, dir));
			assert(seg_id >= 0 && (unsigned)seg_id < segs.size());
			segs.push_back(segs[seg_id]); // clone the segment first
			road_seg_t &seg(segs[seg_id]);
			assert(seg.road_ix < roads.size() && roads[seg.road_ix].dim != dim); // sanity check
			seg        .d[!dim][1] = c.d[!dim][0]; // low part
			segs.back().d[!dim][0] = c.d[!dim][1]; // high part
			cube_t ibc(seg); // intersection bcube
			ibc.d[!dim][0] = c.d[!dim][0]; // copy width from c
			ibc.d[!dim][1] = c.d[!dim][1];
			uint8_t const conns[4] = {7, 11, 13, 14};
			int const other_rix(encode_neg_ix(grn_rix)); // make negative
			isecs[1].emplace_back(ibc, (dim ? seg.road_ix : (int)other_rix), (dim ? other_rix : (int)seg.road_ix), conns[2*(!dim) + dir]);
		}
		float create_connector_road(cube_t const &bcube1, cube_t const &bcube2, vector<cube_t> &blockers, road_network_t *rn1, road_network_t *rn2, unsigned city1, unsigned city2,
			heightmap_query_t &hq, float road_width, float conn_pos, bool dim, bool check_only=0)
		{
			bool const dir(bcube1.d[dim][0] < bcube2.d[dim][0]);
			if (dir == 0) {swap(city1, city2);} // make {lo, hi}
			point p1, p2;
			p1.z = bcube1.d[2][1];
			p2.z = bcube2.d[2][1];
			p1[!dim] = p2[!dim] = conn_pos;
			p1[ dim] = bcube1.d[dim][ dir];
			p2[ dim] = bcube2.d[dim][!dir];
			bool const slope((p1.z < p2.z) ^ dir);
			road_t const road(p1, p2, road_width, dim, slope, roads.size());
			float const road_len(road.get_length()), delta_z(road.get_dz()), max_slope(city_params.max_road_slope);
			assert(road_len > 0.0 && delta_z >= 0.0);
			if (delta_z/road_len > max_slope) {assert(check_only); return -1.0;} // slope is too high (split segments will have even higher slopes)
			unsigned const x1(hq.get_x_pos(road.x1())), y1(hq.get_y_pos(road.y1())), x2(hq.get_x_pos(road.x2())), y2(hq.get_y_pos(road.y2()));

			if (check_only) { // only need to do these checks in this case
				// use find_conn_intersection(4)/make_4way_int() to create a new 4-way intersection on one city if one of these fails?
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
				unsigned const grn_rix(roads.size()); // FIXME: may be wrong end of connector
				if (rn1) {rn1->insert_conn_intersection(road, dim,  dir, grn_rix);}
				if (rn2) {rn2->insert_conn_intersection(road, dim, !dir, grn_rix);}
				float const blocker_padding(max(city_params.road_spacing, 2.0f*city_params.road_border*max(DX_VAL, DY_VAL))); // use road_spacing?
				blockers.push_back(road);
				blockers.back().expand_by(blocker_padding); // add extra padding
			}
			if (road_len <= city_params.conn_road_seg_len) { // simple single road segment case
				if (!check_only) {
					roads.push_back(road);
					road_to_city.emplace_back(city1, city2);
				}
				return hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope]-ROAD_HEIGHT, road.d[2][!slope]-ROAD_HEIGHT, dim, city_params.road_border, check_only);
			}
			unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
			assert(num_segs > 0 && num_segs < 1000); // sanity check
			float const seg_len(road_len/num_segs);
			assert(seg_len <= city_params.conn_road_seg_len);
			road_t rs(road); // keep d[!dim][0], d[!dim][1], dim, and road_ix
			rs.z1() = road.d[2][slope];
			segments.clear();
			float tot_dz(0.0);

			for (unsigned n = 0; n < num_segs; ++n) {
				rs.d[dim][1] = ((n+1 == num_segs) ? road.d[dim][1] : (rs.d[dim][0] + seg_len)); // make sure it ends exactly at the correct location
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
				
				if (!check_only) {
					roads.push_back(*s);
					road_to_city.emplace_back(city1, city2); // Note: city index is specified even for internal (non-terminal) roads
				}
			}
			if (!check_only) { // post-flatten pass to fix up dirt at road joints - doesn't help much
				for (auto s = segments.begin(); s != segments.end(); ++s) {hq.flatten_for_road(*s, city_params.road_border, 0, 1);} // decrease_only=1
			}
			return tot_dz; // success
		}
		void create_connector_bend(cube_t const &int_bcube, bool dx, bool dy, unsigned road_ix_x, unsigned road_ix_y) {
			uint8_t const conns[4] = {6, 5, 10, 9};
			isecs[0].emplace_back(int_bcube, road_ix_x, road_ix_y, conns[2*dy + dx]);
			//blockers.push_back(int_bcube); // ???
		}
		void split_connector_roads(float road_spacing) {
			// Note: here we use segs, maybe 2-way isecs for bends, but not plots
			for (auto r = roads.begin(); r != roads.end(); ++r) {
				bool const d(r->dim), slope(r->slope);
				float const len(r->get_length());
				unsigned const rix(r->road_ix); // not (r - roads.begin())
				if (len <= road_spacing) {segs.emplace_back(*r, rix); continue;} // single segment road
				assert(len > 0.0);
				unsigned const num_segs(ceil(len/road_spacing));
				float const seg_len(len/num_segs), z1(r->d[2][slope]), z2(r->d[2][!slope]); // use fixed-length segments
				assert(seg_len <= road_spacing);
				cube_t c(*r); // start by copying the road's bcube
				
				for (unsigned n = 0; n < num_segs; ++n) {
					c.d[d][1] = ((n+1 == num_segs) ? r->d[d][1] : (c.d[d][0] + seg_len)); // make sure it ends exactly at the correct location
					for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
					if (c.z2() < c.z1()) {swap(c.z2(), c.z1());} // swap zvals if needed
					assert(c.is_normalized());
					segs.emplace_back(c, rix, d, r->slope);
					c.d[d][0] = c.d[d][1]; // shift segment end point
				} // for n
			} // for r
		}
		void gen_tile_blocks(road_network_t const &global_rn) {
			tile_blocks.clear(); // should already be empty?
			map<uint64_t, unsigned> tile_to_block_map;
			add_tile_blocks(segs,  tile_to_block_map, TYPE_RSEG);
			add_tile_blocks(plots, tile_to_block_map, TYPE_PLOT);
			for (unsigned i = 0; i < 3; ++i) {add_tile_blocks(isecs[i], tile_to_block_map, (TYPE_ISEC2 + i));}
			//cout << "tile_to_block_map: " << tile_to_block_map.size() << ", tile_blocks: " << tile_blocks.size() << endl;
			calc_ix_values(global_rn);
		}
		void get_road_bcubes(vector<cube_t> &bcubes) const {get_all_bcubes(roads, bcubes);}
		void get_plot_bcubes(vector<cube_t> &bcubes) const {get_all_bcubes(plots, bcubes);}

		bool check_road_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only=1) const {
			if (roads.empty()) return 0;
			point const query_pos(pos - get_query_xlate());
			if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
			if (check_bcubes_sphere_coll(roads, query_pos, radius, xy_only)) return 1;

			if (include_intersections) { // used for global road network
				for (unsigned i = 0; i < 3; ++i) { // {2-way, 3-way, 4-way}
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
		static float get_car_lane_offset() {return CAR_LANE_OFFSET*city_params.road_width;}

		bool add_car(car_t &car, rand_gen_t &rgen) const {
			if (segs.empty()) return 0; // no segments to place car on

			for (unsigned n = 0; n < 10; ++n) { // make 10 tries
				unsigned const seg_ix(rgen.rand()%segs.size());
				road_seg_t const &seg(segs[seg_ix]); // chose a random segment
				car.dim   = seg.dim;
				car.dir   = (rgen.rand()&1);
				car.dz    = 0.0; // flat
				car.max_speed = rgen.rand_uniform(0.66, 1.0); // add some speed variation
				car.cur_speed = 0.0;
				car.cur_road  = seg.road_ix;
				car.cur_seg   = seg_ix;
				car.cur_road_type = TYPE_RSEG;
				car.turn_dir = TURN_NONE;
				vector3d car_sz; // {length, width, height}
				for (unsigned d = 0; d < 3; ++d) {car_sz[d] = CAR_SIZE[d]*rgen.rand_uniform(0.9, 1.1)*city_params.road_width;}
				car.height = car_sz.z;
				point pos;
				float val1(seg.d[seg.dim][0] + 0.5*car_sz.x), val2(seg.d[seg.dim][1] - 0.5*car_sz.x);
				if (val1 >= val2) continue; // failed, try again (connector road junction?)
				pos[!seg.dim]  = 0.5*(seg.d[!seg.dim][0] + seg.d[!seg.dim][1]); // center of road
				pos[!seg.dim] += ((car.dir ^ car.dim) ? -1.0 : 1.0)*get_car_lane_offset(); // place in right lane
				pos[ seg.dim] = rgen.rand_uniform(val1, val2); // place at random pos in segment
				pos.z = seg.z2() + 0.5*car_sz.z; // place above road surface
				if (seg.dim) {swap(car_sz.x, car_sz.y);}
				car.bcube.set_from_point(pos);
				car.bcube.expand_by(0.5*car_sz);
				assert(get_road_bcube_for_car(car).contains_cube_xy(car.bcube)); // sanity check
				return 1; // success
			} // for n
			return 0; // failed
		}
		bool find_car_next_seg(car_t &car, road_network_t const &global_rn) const {
			if (car.cur_road_type == TYPE_RSEG) {
				road_seg_t const &seg(get_car_seg(car));
				//cout << "before: " << car.str() << endl << "before road bcube: " << get_road_bcube_for_car(car).str() << endl;
				car.cur_road_type = seg.conn_type[car.dir];
				car.cur_road      = seg.road_ix;
				car.cur_seg       = seg.conn_ix[car.dir];
				//car.cur_city = ?; // FIXME: transition from global road network back to local city roads
				//cout << "after : " << car.str() << endl << "after  road bcube: " << get_road_bcube_for_car(car).str() << endl;
				assert(get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return 1; // always within same city, no city update
			}
			road_isec_t const &isec(get_car_isec(car)); // conn_ix: {-x, +x, -y, +y}
			unsigned const orient(car.get_orient());
			assert(isec.conn & (1<<orient));
			int conn_ix(isec.conn_ix[orient]), rix(isec.rix_xy[car.dim]);

			if (conn_ix < 0) { // city connector road case, use global_rn
				assert(rix < 0);
				conn_ix = decode_neg_ix(conn_ix);
				assert((unsigned)conn_ix < global_rn.segs.size());
				rix = global_rn.segs[conn_ix].road_ix;
				car.cur_city = CONN_CITY_IX; // move to global road network
			}
			else {
				assert(rix >= 0); // local road
			}
			car.cur_road = (unsigned)rix;
			car.cur_seg  = (unsigned)conn_ix;
			car.cur_road_type = TYPE_RSEG; // always connects to a road segment
			cube_t const bcube(get_road_bcube_for_car(car, global_rn));

			if (!bcube.intersects_xy(car.bcube)) { // sanity check
				cout << "bad intersection:" << endl << car.str() << endl << bcube.str() << endl;
				assert(0);
			}
			return 1; // done
		}
		void update_car(car_t &car, rand_gen_t &rgen, road_network_t const &global_rn) const {
			if (car.is_parked()) return; // stopped, no update (for now)

			if (car.stopped_at_light) { // Note: is_isect test is here to allow cars to coast through lights when decel is very low
				bool const was_stopped(car.cur_speed == 0.0);
				if (!is_isect(car.cur_road_type) || !get_car_isec(car).red_light(car)) {car.stopped_at_light = 0;} // can go now
				else {car.decelerate_fast();}
				if (was_stopped) return; // no update needed
			} else {car.accelerate();}
			cube_t const bcube(get_road_bcube_for_car(car));
			if (!bcube.intersects_xy(car.prev_bcube)) {cout << car.str() << endl << bcube.str() << endl; assert(0);} // sanity check
			unsigned conn_left[4] = {3,2,0,1}, conn_right[4] = {2,3,1,0};
			bool const dim(car.dim);
			float const road_dz(bcube.get_dz());

			if (car.cur_road_type == TYPE_RSEG && road_dz != 0.0) {
				assert(car.cur_city == CONN_CITY_IX);
				bool const slope(get_car_seg(car).slope);
				float const car_pos(0.5*(car.bcube.d[dim][0] + car.bcube.d[dim][1])); // center of car in dim
				float const road_len(bcube.d[dim][1] - bcube.d[dim][0]);
				float const t((car_pos - bcube.d[dim][0])/road_len); // car pos along road in (0.0, 1.0)
				float const road_z(bcube.d[2][slope] + t*(bcube.d[2][!slope] - bcube.d[2][slope]));
				float const car_len(car.get_length());
				car.dz = ((slope ^ car.dir) ? 1.0 : -1.0)*road_dz*(car_len/road_len);
				car.bcube.z1() = road_z - 0.5*fabs(car.dz);
				car.bcube.z2() = road_z + 0.5*fabs(car.dz) + car.height;
			}
			if (bcube.contains_cube_xy(car.bcube)) { // in same road seg/int
				if (car.turn_dir != TURN_NONE) {
					assert(is_isect(car.cur_road_type));
					float const car_lane_offset(get_car_lane_offset()), isec_center(bcube.get_cube_center()[dim]);
					float const centerline(isec_center + (((car.turn_dir == TURN_LEFT) ^ car.dir) ? -1.0 : 1.0)*car_lane_offset);
					float const prev_val(car.prev_bcube.get_cube_center()[dim]), cur_val(car.bcube.get_cube_center()[dim]);
					
					if (min(prev_val, cur_val) <= centerline && max(prev_val, cur_val) > centerline) { // crossed the lane centerline boundary
						car.move_by(centerline - cur_val); // align to lane centerline
						vector3d const car_sz(car.bcube.get_size());
						float const size_adj(0.5*(car_sz[dim] - car_sz[!dim]));
						vector3d expand(zero_vector);
						expand[dim] -= size_adj; expand[!dim] += size_adj;
						car.bcube.expand_by(expand); // fix aspect ratio
						if ((dim == 0) ^ (car.turn_dir == TURN_LEFT)) {car.dir ^= 1;}
						car.dim ^= 1;
						car.turn_dir = TURN_NONE; // turn completed
					}
				}
				assert(get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return; // done
			}
			// car crossing the border of this bcube, update state
			if (!bcube.contains_pt_xy_inc_low_edge(car.bcube.get_cube_center())) { // move to another road seg/int
				if (!find_car_next_seg(car, global_rn)) {car.park(); return;} // update failed, park/stop car (for now)

				if (is_isect(car.cur_road_type)) { // moved into an intersection, choose direction
					road_isec_t const &isec(get_car_isec(car));
					unsigned const orient_in(2*car.dim + (!car.dir)); // invert dir (incoming, not outgoing)
					assert(isec.conn & (1<<orient_in)); // car must come from an enabled orient
					unsigned orients[3]; // {straight, left, right}
					orients[TURN_NONE ] = car.get_orient(); // straight
					orients[TURN_LEFT ] = conn_left [orient_in];
					orients[TURN_RIGHT] = conn_right[orient_in];
					
					// TODO: path finding update - use A*?
					while (1) {
						unsigned new_turn_dir(0);
						int const rval(rand()%4);
						if      (rval == 0) {new_turn_dir = TURN_LEFT ;} // 25%
						else if (rval == 1) {new_turn_dir = TURN_RIGHT;} // 25%
						else                {new_turn_dir = TURN_NONE ;} // 50%
						if (isec.conn & (1<<orients[new_turn_dir])) {car.turn_dir = new_turn_dir; break;}
					} // end while
					//cout << TXT(orient_in) << TXT(car.get_orient()) << "turn_dir=" << unsigned(car.turn_dir) << " conn=" << unsigned(isec.conn) << endl;
					assert(isec.conn & (1<<orients[car.turn_dir]));
					car.stopped_at_light = isec.red_or_yellow_light(car); // FIXME: check this earlier, before the car is in the intersection?
					if (car.stopped_at_light) {car.decelerate_fast();}
				}
			}
			assert(get_road_bcube_for_car(car, global_rn).intersects_xy(car.bcube)); // sanity check
		}
		void next_frame() {
			for (unsigned n = 1; n < 3; ++n) { // {2-way, 3-way, 4-way} - Note: 2-way can be skipped
				for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {i->next_frame();} // update stoplight state
			}
		}
		road_seg_t const &get_car_seg(car_t const &car) const {
			assert(car.cur_road_type == TYPE_RSEG);
			assert(car.cur_seg < segs.size());
			return segs[car.cur_seg];
		}
		road_isec_t const &get_car_isec(car_t const &car) const {
			assert(is_isect(car.cur_road_type));
			auto const &iv(isecs[car.cur_road_type - TYPE_ISEC2]);
			assert(car.cur_seg < iv.size());
			return iv[car.cur_seg];
		}
		cube_t get_road_bcube_for_car(car_t const &car) const {
			assert(car.cur_road < roads.size());
			if (car.cur_road_type == TYPE_RSEG) {return get_car_seg(car);}
			return get_car_isec(car);
		}
		cube_t get_road_bcube_for_car(car_t const &car, road_network_t const &global_rn) const { // function variant that works with both global and local roads
			if (car.cur_city == city_id) {return get_road_bcube_for_car(car);}
			else if (car.cur_city == CONN_CITY_IX) {return global_rn.get_road_bcube_for_car(car);}
			else {assert(0); return cube_t();}
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
		road_networks.push_back(road_network_t(region, road_networks.size()));
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
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, hq, road_width, conn_pos, !d, 1)); // check_only=1
					if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost;}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Single segment dim: << "d " << cost: " << best_cost << endl;
					global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, hq, road_width, best_conn_pos, !d, 0); // check_only=0; actually make the change
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
					float const cost1(global_rn.create_connector_road(bcube1, int_cube, blockers, &rn1, nullptr, city1, CONN_CITY_IX, hq, road_width, (fdim ? xval : yval),  fdim, 1)); // check_only=1
					if (cost1 < 0.0) continue; // bad segment
					if (best_cost > 0.0 && cost1 > best_cost) continue; // bound - early terminate
					float const cost2(global_rn.create_connector_road(int_cube, bcube2, blockers, nullptr, &rn2, CONN_CITY_IX, city2, hq, road_width, (fdim ? yval : xval), !fdim, 1)); // check_only=1
					if (cost2 < 0.0) continue; // bad segment
					float const cost(cost1 + cost2); // Note: cost function will prefer shorter routes
					if (best_cost < 0.0 || cost < best_cost) {best_xval = xval; best_yval = yval; best_int_cube = int_cube; best_cost = cost;}
				} // for n
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Double segment cost: " << best_cost << " " << TXT(best_xval) << TXT(best_yval) << TXT(fdim) << ", int_cube: " << best_int_cube.str() << endl;
					hq.flatten_region_to(best_int_cube, city_params.road_border); // do this first to improve flattening
					unsigned road_ix[2];
					road_ix[ fdim] = global_rn.num_roads();
					global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, city1, CONN_CITY_IX, hq, road_width, (fdim ? best_xval : best_yval),  fdim, 0); // check_only=0
					flatten_op_t const fop(hq.last_flatten_op); // cache for reuse later during decrease_only pass
					road_ix[!fdim] = global_rn.num_roads();
					global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2, CONN_CITY_IX, city2, hq, road_width, (fdim ? best_yval : best_xval), !fdim, 0); // check_only=0
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), road_ix[0], road_ix[1]);
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
		timer_t timer("Gen Tile Blocks");
		global_rn.gen_tile_blocks(global_rn); // must be done first to fill in road_to_city and city_to_seg
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_tile_blocks(global_rn);}
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
	void next_frame() {
		//timer_t timer("Update Stoplights");
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->next_frame();}
		global_rn.next_frame(); // not needed since there are no 3/4-way intersections/stoplights?
	}
	bool add_car(car_t &car, rand_gen_t &rgen) const {
		if (road_networks.empty()) return 0; // no cities to add cars to
		unsigned const city(rgen.rand()%road_networks.size());
		car.cur_city = city;
		return road_networks[city].add_car(car, rgen);
	}
	road_network_t const &get_car_rn(car_t const &car) const {
		if (car.cur_city == CONN_CITY_IX) {return global_rn;}
		assert(car.cur_city < road_networks.size());
		return road_networks[car.cur_city];
	}
	void update_car(car_t &car, rand_gen_t &rgen) const {get_car_rn(car).update_car(car, rgen, global_rn);}
	cube_t get_road_bcube_for_car(car_t const &car) const {return get_car_rn(car).get_road_bcube_for_car(car);}
}; // city_road_gen_t


bool car_t::check_collision(car_t &c, city_road_gen_t const &road_gen) {

	//if (dim != c.dim || dir != c.dir) return 0;
	float const avg_len(0.5*((bcube.d[dim][1] - bcube.d[dim][0]) + (c.bcube.d[c.dim][1] - c.bcube.d[c.dim][0]))); // average length of the two cars
	float const min_speed(min(cur_speed, c.cur_speed)); // relative to max speed of 1.0
	float const sep_dist(avg_len*(0.25 + 1.0*min_speed)); // 25% to 125% car length, depending on speed
	float const test_dist(0.999*sep_dist); // slightly smaller than separation distance
	cube_t bcube_ext(bcube);
	bcube_ext.d[dim][0] -= test_dist; bcube_ext.d[dim][1] += test_dist; // expand by test_dist distance
	if (!bcube_ext.intersects(c.bcube)) return 0;
	if (c.dim != dim) return 0; // turning in an intersection
	assert(c.dir == dir); // otherwise should be on different sides of the road
	float const front(bcube.d[dim][dir]), c_front(c.bcube.d[c.dim][c.dir]);
	bool const move_c((front < c_front) ^ dir); // move the car that's behind
	// Note: we could slow the car in behind, but that won't work for initial placement collisions when speed == 0
	car_t &cmove(move_c ? c : *this); // the car that will be moved
	car_t const &cstay(move_c ? *this : c); // the car that won't be moved
	//cout << "Collision between " << cmove.str() << " and " << cstay.str() << endl;
	if (cstay.is_stopped()) {cmove.decelerate_fast();} else {cmove.decelerate();}
	float const dist(cstay.bcube.d[dim][!dir] - cmove.bcube.d[dim][dir]); // signed distance between the back of the car in front, and the front of the car in back
	point delta(all_zeros);
	delta[dim] += dist + (cmove.dir ? -sep_dist : sep_dist); // force separation between cars
	cube_t const &bcube(road_gen.get_road_bcube_for_car(cmove));

	if (!bcube.contains_cube_xy(cmove.bcube + delta)) { // moved outside its current road segment bcube
		//if (cmove.bcube == cmove.prev_bcube) {return 1;} // collided, but not safe to move the car (init pos or second collision)
		if (cmove.bcube != cmove.prev_bcube) { // try resetting to last frame's position
			cmove.bcube  = cmove.prev_bcube; // restore prev frame's pos
			return 1; // done
		}
		else { // keep the car from moving outside its current segment (init collision case)
			if (cmove.dir) {max_eq(delta[dim], min(0.0f, 0.999f*(bcube.d[cmove.dim][0] - cmove.bcube.d[cmove.dim][0])));}
			else           {min_eq(delta[dim], max(0.0f, 0.999f*(bcube.d[cmove.dim][1] - cmove.bcube.d[cmove.dim][1])));}
		}
	}
	cmove.bcube += delta;
	return 1;
}


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

		void draw_car(car_t const &car) { // Note: all quads
			if (!check_cube_visible(car.bcube)) return;
			begin_tile(car.get_center()); // enable shadows
			quad_batch_draw &qbd(qbds[emit_now]);
			assert(car.color_id < NUM_CAR_COLORS);
			colorRGBA const &color(car_colors[car.color_id]);
			cube_t const &c(car.bcube);
			bool const d(car.dim), D(car.dir);
			float const dz(car.dz);//, dx(dz*car.height/car.get_length());
			float fz1, fz2, bz1, bz2; // front/back top/bottom
			if (dz == 0.0) {} // FIXME: optimization for level car case (within city)?
			if (dz > 0.0) {fz1 = c.z1() + dz; fz2 = c.z2(); bz1 = c.z1(); bz2 = c.z2() - dz;} // going uphill
			else          {fz1 = c.z1(); fz2 = c.z2() + dz; bz1 = c.z1() - dz; bz2 = c.z2();} // going downhill or level
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
			if (emit_now) {qbds[1].draw_and_clear();} // shadowed (FIXME: only when tile changes)
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
		if (num == 0) return;
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
		if (cars.empty() || !animate2) return;
		//timer_t timer("Update Cars");
		sort(cars.begin(), cars.end()); // sort by city/road/segment; is this a good idea?
		float const speed(0.001*car_speed*fticks);
		for (auto i = cars.begin(); i != cars.end(); ++i) {i->move(speed);} // move cars

		for (auto i = cars.begin(); i != cars.end(); ++i) { // collision detection
			for (auto j = i+1; j != cars.end(); ++j) { // check for collisions with cars on the same road (can't test seg because they can be on diff segs but still collide)
				if (i->cur_city == j->cur_city && i->cur_road == j->cur_road) {i->check_collision(*j, road_gen);}
				else {break;}
			}
		} // for i
		for (auto i = cars.begin(); i != cars.end(); ++i) {road_gen.update_car(*i, rgen);} // run update logic
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
	void next_frame() {
		road_gen.next_frame(); // update stoplights
		car_manager.next_frame(city_params.car_speed);
	}
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

