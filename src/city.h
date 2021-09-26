// 3D World - City Header
// by Frank Gennari
// 11/19/18

#pragma once

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"
#include "shaders.h"
#include "draw_utils.h"
#include "buildings.h" // for building_occlusion_state_t and obj models
#include "city_model.h"

using std::string;

unsigned const CONN_CITY_IX((1<<16)-1); // uint16_t max
unsigned const NO_CITY_IX(CONN_CITY_IX-1); // used for cars not in any city (in house garages)

enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   TID_PARK_LOT,  TID_TRACKS,  TID_PARK,  TID_DRIVEWAY,  /*none for bldg*/ NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, TYPE_PARK_LOT, TYPE_TRACKS, TYPE_PARK, TYPE_DRIVEWAY, TYPE_BUILDING,    NUM_RD_TYPES};
enum {TURN_NONE=0, TURN_LEFT, TURN_RIGHT, TURN_UNSPEC};
enum {INT_NONE=0, INT_ROAD, INT_PLOT, INT_PARKING, INT_PARK};
enum {RTYPE_ROAD=0, RTYPE_TRACKS};
unsigned const CONN_TYPE_NONE = 0;
colorRGBA const road_colors[NUM_RD_TYPES] = {WHITE, WHITE, WHITE, WHITE, WHITE, WHITE, WHITE, LT_GRAY, LT_GRAY, WHITE}; // all white except for parks and driveways

int       const FORCE_MODEL_ID = -1; // -1 disables
unsigned  const NUM_CAR_COLORS = 10;
colorRGBA const car_colors[NUM_CAR_COLORS] = {WHITE, GRAY_BLACK, GRAY, ORANGE, RED, DK_RED, DK_BLUE, colorRGBA(0.5, 0.9, 0.5), YELLOW, BROWN};

float const ROAD_HEIGHT          = 0.002;
float const PARK_SPACE_WIDTH     = 1.6;
float const PARK_SPACE_LENGTH    = 1.8;
float const CONN_ROAD_SPEED_MULT = 2.0; // twice the speed limit on connector roads
float const HEADLIGHT_ON_RAND    = 0.1;
float const STREETLIGHT_ON_RAND  = 0.05;
float const TUNNEL_WALL_THICK    = 0.25; // relative to radius
float const TRACKS_WIDTH         = 0.5; // relative to road width
float const CAR_SPEED_SCALE      = 0.001;
float const SIDEWALK_WIDTH       = 0.1; // relative to road texture
vector3d const CAR_SIZE(0.30, 0.13, 0.08); // {length, width, height} in units of road width
float const CAR_RADIUS_SCALE(CAR_SIZE.mag()/CAR_SIZE.z);

float const PED_WIDTH_SCALE  = 0.5; // ratio of collision radius to model radius (x/y)
float const PED_HEIGHT_SCALE = 2.5; // ratio of collision radius to model height (z)

extern double tfticks;
extern float fticks;


inline bool is_isect(unsigned type) {return (type >= TYPE_ISEC2 && type <= TYPE_ISEC4);}
inline int encode_neg_ix(unsigned ix) {return -(int(ix)+1);}
inline unsigned decode_neg_ix(int ix) {assert(ix < 0); return -(ix+1);}

inline float rand_hash(float to_hash) {return fract(12345.6789*to_hash);}
inline float signed_rand_hash(float to_hash) {return 0.5*(rand_hash(to_hash) - 1.0);}


struct city_params_t {

	unsigned num_cities, num_samples, num_conn_tries, city_size_min, city_size_max, city_border, road_border, slope_width, num_rr_tracks, park_rate;
	float road_width, road_spacing, road_spacing_rand, road_spacing_xy_add, conn_road_seg_len, max_road_slope, residential_probability;
	unsigned make_4_way_ints; // 0=all 3-way intersections; 1=allow 4-way; 2=all connector roads must have at least a 4-way on one end; 3=only 4-way (no straight roads)
	bool assign_house_plots, new_city_conn_road_alg;
	// cars
	unsigned num_cars;
	float car_speed, traffic_balance_val, new_city_prob, max_car_scale;
	bool enable_car_path_finding, convert_model_files;
	vector<city_model_t> car_model_files, ped_model_files, hc_model_files;
	city_model_t fire_hydrant_model_file;
	// parking lots
	unsigned min_park_spaces, min_park_rows;
	float min_park_density, max_park_density;
	// lighting
	bool car_shadows;
	unsigned max_lights, max_shadow_maps, smap_size;
	// trees
	unsigned max_trees_per_plot;
	float tree_spacing;
	// detail objects
	unsigned max_benches_per_plot;
	// pedestrians
	unsigned num_peds, num_building_peds;
	float ped_speed;
	bool ped_respawn_at_dest;
	// buildings; maybe should be building params, but we have the model loading code here
	city_model_t building_models[NUM_OBJ_MODELS];

	city_params_t() : num_cities(0), num_samples(100), num_conn_tries(50), city_size_min(0), city_size_max(0), city_border(0), road_border(0), slope_width(0),
		num_rr_tracks(0), park_rate(0), road_width(0.0), road_spacing(0.0), road_spacing_rand(0.0), road_spacing_xy_add(0.0), conn_road_seg_len(1000.0),
		max_road_slope(1.0), residential_probability(0.0), make_4_way_ints(0), assign_house_plots(0), new_city_conn_road_alg(0), num_cars(0), car_speed(0.0),
		traffic_balance_val(0.5), new_city_prob(1.0), max_car_scale(1.0), enable_car_path_finding(0), convert_model_files(0), min_park_spaces(12),
		min_park_rows(1), min_park_density(0.0), max_park_density(1.0), car_shadows(0), max_lights(1024), max_shadow_maps(0), smap_size(0),
		max_trees_per_plot(0), tree_spacing(1.0), max_benches_per_plot(0), num_peds(0), num_building_peds(0), ped_speed(0.0), ped_respawn_at_dest(0) {}
	bool enabled() const {return (num_cities > 0 && city_size_min > 0);}
	bool roads_enabled() const {return (road_width > 0.0 && road_spacing > 0.0);}
	float get_road_ar() const {return round(road_spacing/road_width);} // round to nearest texture multiple
	static bool read_error(string const &str) {cout << "Error reading city config option " << str << "." << endl; return 0;}
	bool read_option(FILE *fp);
	bool add_model(unsigned id, FILE *fp);
	bool has_helicopter_model() const {return !hc_model_files.empty();}
	vector3d get_nom_car_size() const {return CAR_SIZE*road_width;}
	vector3d get_max_car_size() const {return max_car_scale*get_nom_car_size();}
}; // city_params_t


struct car_base_t;

struct road_gen_base_t {
	virtual cube_t get_bcube_for_car(car_base_t const &car) const = 0;
 	virtual ~road_gen_base_t() {}
};


struct waiting_obj_t {
	float waiting_start;
	waiting_obj_t() : waiting_start(0.0) {}
	void reset_waiting() {waiting_start = tfticks;}
	float get_wait_time_secs() const {return (float(tfticks) - waiting_start)/TICKS_PER_SECOND;} // Note: only meaningful for cars stopped at lights or peds stopped at roads
};

struct car_base_t { // the part needed for the pedestrian interface (size = 36)
	cube_t bcube;
	bool dim, dir, stopped_at_light;
	unsigned char cur_road_type, turn_dir;
	unsigned short cur_city, cur_road, cur_seg;
	float max_speed, cur_speed;

	car_base_t() : bcube(all_zeros), dim(0), dir(0), stopped_at_light(0), cur_road_type(TYPE_RSEG), turn_dir(TURN_NONE), cur_city(0), cur_road(0), cur_seg(0), max_speed(0.0), cur_speed(0.0) {}
	point get_center() const {return bcube.get_cube_center();}
	unsigned get_orient() const {return (2*dim + dir);}
	unsigned get_orient_in_isec() const {return (2*dim + (!dir));} // invert dir (incoming, not outgoing)
	bool in_isect() const {return is_isect(cur_road_type);}
	unsigned get_isec_type() const {assert(in_isect()); return (cur_road_type - TYPE_ISEC2);}
	float get_max_speed() const {return ((cur_city == CONN_CITY_IX) ? CONN_ROAD_SPEED_MULT : 1.0)*max_speed;}
	float get_length() const {return bcube.get_sz_dim( dim);}
	float get_width () const {return bcube.get_sz_dim(!dim);}
	bool is_almost_stopped() const {return (cur_speed < 0.1*max_speed);}
	bool is_stopped () const {return (cur_speed == 0.0);}
	bool is_parked  () const {return (max_speed == 0.0);}
	bool in_garage  () const {return (cur_city == NO_CITY_IX);}
	point get_front(float dval=0.5) const;
};

struct car_t : public car_base_t, public waiting_obj_t { // size = 96
	cube_t prev_bcube;
	bool entering_city, in_tunnel, dest_valid, destroyed, in_reverse;
	unsigned char color_id, front_car_turn_dir, model_id;
	unsigned short dest_city, dest_isec;
	short dest_driveway; // -1 is unset
	float height, dz, rot_z, turn_val, waiting_pos;
	car_t const *car_in_front;

	car_t() : prev_bcube(all_zeros), entering_city(0), in_tunnel(0), dest_valid(0), destroyed(0), in_reverse(0), color_id(0), front_car_turn_dir(TURN_UNSPEC),
		model_id(0), dest_city(0), dest_isec(0), dest_driveway(-1), height(0.0), dz(0.0), rot_z(0.0), turn_val(0.0), waiting_pos(0.0), car_in_front(nullptr) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	float get_max_lookahead_dist() const;
	bool headlights_on() const;
	float get_turn_rot_z(float dist_to_turn) const;
	colorRGBA const &get_color() const {assert(color_id < NUM_CAR_COLORS); return car_colors[color_id];}
	void apply_scale(float scale);
	void destroy();
	float get_min_sep_dist_to_car(car_t const &c, bool add_one_car_len=0) const;
	string str() const;
	string label_str() const;
	void move(float speed_mult);
	void maybe_accelerate(float mult=0.02);
	void accelerate(float mult=0.02) {cur_speed = min(get_max_speed(), (cur_speed + mult*fticks*max_speed));}
	void decelerate(float mult=0.05) {cur_speed = max(0.0f, (cur_speed - mult*fticks*max_speed));}
	void decelerate_fast() {decelerate(10.0);} // Note: large decel to avoid stopping in an intersection
	void park() {cur_speed = max_speed = 0.0;}
	void stop() {cur_speed = 0.0;} // immediate stop
	void move_by(float val) {bcube.translate_dim(dim, val);}
	void begin_turn() {turn_val = bcube.get_center_dim(!dim);}
	bool maybe_apply_turn(float centerline, bool for_driveway);
	bool check_collision(car_t &c, road_gen_base_t const &road_gen);
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm) const;
	bool front_intersects_car(car_t const &c) const;
	void honk_horn_if_close() const;
	void honk_horn_if_close_and_fast() const;
	void on_alternate_turn_dir(rand_gen_t &rgen);
	void register_adj_car(car_t &c);
	unsigned count_cars_in_front(cube_t const &range=cube_t(all_zeros)) const;
	float get_sum_len_space_for_cars_in_front(cube_t const &range) const;
};

struct car_city_vect_t {
	vector<car_base_t> cars[2][2]; // {dim x dir}
	vector<cube_with_ix_t> parked_car_bcubes; // stores car bcube + plot_ix
	void clear_cars();
	void clear() {clear_cars(); parked_car_bcubes.clear();}
};


struct comp_car_road {
	bool operator()(car_base_t const &c1, car_base_t const &c2) const {return (c1.cur_road < c2.cur_road);}
};
struct comp_car_city_then_road {
	bool operator()(car_base_t const &c1, car_base_t const &c2) const {
		return ((c1.cur_city != c2.cur_city) ? (c1.cur_city < c2.cur_city) : (c1.cur_road < c2.cur_road));
	}
};
struct comp_car_road_then_pos {
	vector3d const &camera_pos;
	comp_car_road_then_pos(vector3d const &camera_pos_) : camera_pos(camera_pos_) {}
	bool operator()(car_t const &c1, car_t const &c2) const;
};


struct helicopter_t {
	enum {STATE_WAIT=0, STATE_TAKEOFF, STATE_FLY, STATE_LAND};
	cube_t bcube;
	vector3d dir;
	// dynamic state for moving helicopters
	vector3d velocity;
	float wait_time; // time to wait before takeoff
	float fly_zval; // zval required for flight to avoid buildings and terrain
	float blade_rot;
	unsigned dest_hp; // destination (or current) helipad
	unsigned state, model_id;
	bool dynamic, dynamic_shadow;

	helicopter_t(cube_t const &bcube_, vector3d const &dir_, unsigned model_id_, unsigned dest_hp_, bool dynamic_) :
		bcube(bcube_), dir(dir_), velocity(zero_vector), wait_time(0.0), fly_zval(0.0), blade_rot(0.0),
		dest_hp(dest_hp_), state(STATE_WAIT), model_id(model_id_), dynamic(dynamic_), dynamic_shadow(0) {}
	point get_landing_pt() const {return point(bcube.xc(), bcube.yc(), bcube.z1());}
	void invalidate_tile_shadow_map(vector3d const &shadow_offset, bool repeat_next_frame) const;
};


class road_mat_mgr_t {

	bool inited;
	unsigned tids[NUM_RD_TIDS] = {}, sl_tid;

public:
	road_mat_mgr_t() : inited(0), sl_tid(0) {}
	void ensure_road_textures();
	void set_texture(unsigned type);
	void set_stoplight_texture();
};

template<typename T> static void add_flat_city_quad(T const &r, quad_batch_draw &qbd, colorRGBA const &color, float ar) { // z1 == z2
	float const z(r.z1());
	point const pts[4] = {point(r.x1(), r.y1(), z), point(r.x2(), r.y1(), z), point(r.x2(), r.y2(), z), point(r.x1(), r.y2(), z)};
	qbd.add_quad_pts(pts, color, plus_z, r.get_tex_range(ar));
}

struct road_t : public cube_t {
	unsigned road_ix;
	//unsigned char type; // road, railroad, etc. {RTYPE_ROAD, RTYPE_TRACKS}
	bool dim; // dim the road runs in
	bool slope; // 0: z1 applies to first (lower) point; 1: z1 applies to second (upper) point

	road_t(cube_t const &c, bool dim_, bool slope_=0, unsigned road_ix_=0) : cube_t(c), road_ix(road_ix_), dim(dim_), slope(slope_) {}
	road_t(point const &s, point const &e, float width, bool dim_, bool slope_=0, unsigned road_ix_=0);
	float get_length   () const {return (d[ dim][1] - d[ dim][0]);}
	float get_width    () const {return (d[!dim][1] - d[!dim][0]);}
	float get_slope_val() const {return dz()/get_length();}
	float get_start_z  () const {return (slope ? z2() : z1());}
	float get_end_z    () const {return (slope ? z1() : z2());}
	float get_z_adj    () const {return (ROAD_HEIGHT + 0.5*get_slope_val()*(dim ? DY_VAL : DX_VAL));} // account for a half texel of error for sloped roads
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, -ar, (dim ? -1.0 : 1.0), 0, dim);}
	cube_t const &get_bcube() const {return *this;}
	cube_t       &get_bcube()       {return *this;}
	void add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar) const;
};

struct road_seg_t : public road_t {
	unsigned short road_ix, conn_ix[2], conn_type[2]; // {dim=0, dim=1}
	mutable unsigned short car_count; // can be written to during car update logic

	void init_ixs() {conn_ix[0] = conn_ix[1] = 0; conn_type[0] = conn_type[1] = CONN_TYPE_NONE;}
	road_seg_t(road_t const &r, unsigned rix) : road_t(r), road_ix(rix), car_count(0) {init_ixs();}
	road_seg_t(cube_t const &c, unsigned rix, bool dim_, bool slope_=0) : road_t(c, dim_, slope_), road_ix(rix), car_count(0) {init_ixs();}
	void next_frame() {car_count = 0;}
};

struct driveway_t : public cube_t {
	bool dim, dir; // direction to road; d[dim][dir] is the edge shared with the road
	mutable bool in_use; // either reserves the spot, or car <car_ix> is parked there; modified by car_manager in a different thread - must be mutable, maybe should be atomic
	unsigned plot_ix;
	int car_ix; // driveway can old exactly one car; -1 = none
	driveway_t() : dim(0), dir(0), in_use(0), plot_ix(0), car_ix(-1) {}
	driveway_t(cube_t const &c, bool dim_, bool dir_, unsigned pix) : cube_t(c), dim(dim_), dir(dir_), in_use(0), plot_ix(pix), car_ix(-1) {}
	float get_edge_at_road() const {return d[dim][dir];}
	float get_length() const {return get_sz_dim(dim);}
	int get_cur_car()  const {return (in_use ? car_ix : -1);} // if not in use, return -1/none even if car_ix hasn't been reset
	void add_car(unsigned ix) {assert(!in_use); car_ix = ix; in_use = 1;}
	tex_range_t get_tex_range(float ar) const;
};

struct road_plot_t : public cube_t {
	uint8_t xpos, ypos; // position within the city grid
	bool is_residential, has_parking, is_park;
	road_plot_t(cube_t const &c, uint8_t xpos_, uint8_t ypos_, bool is_res=0) : cube_t(c), xpos(xpos_), ypos(ypos_), is_residential(is_res), has_parking(0), is_park(0) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, ar, ar);}
};

struct plot_adj_t {
	int adj[4]; // adjacent plots, by index: west, east, south, north; -1 is none
	plot_adj_t() {adj[0] = adj[1] = adj[2] = adj[3] = -1;}
	void set_adj(unsigned dir, unsigned ix) {assert(dir < 4); assert(adj[dir] == -1); adj[dir] = ix;} // can only set once
};

struct plot_xy_t {
	unsigned nx, ny;
	vector<plot_adj_t> adj_plots;
	plot_xy_t() : nx(0), ny(0) {}
	plot_xy_t(unsigned nx_, unsigned ny_) : nx(nx_), ny(ny_) {}
	void gen_adj_plots(vector<road_plot_t> const &plots);
	unsigned num() const {return nx*ny;}
	int get_adj(unsigned x, unsigned y, unsigned dir) const {assert(x < nx && y < ny); return adj_plots[y*nx + x].adj[dir];} // -1 == no plot
};

struct parking_lot_t : public cube_t {
	bool dim, dir;
	unsigned short row_sz, num_rows;
	vector<unsigned char> used_spaces;
	parking_lot_t(cube_t const &c, bool dim_, bool dir_, unsigned row_sz_=0, unsigned num_rows_=0) : cube_t(c), dim(dim_), dir(dir_), row_sz(row_sz_), num_rows(num_rows_) {}
	tex_range_t get_tex_range(float ar) const;
};

namespace stoplight_ns {

	enum {GREEN_LIGHT=0, YELLOW_LIGHT, RED_LIGHT}; // colors, unused (only have stop and go states anyway)
	enum {EGL=0, EGWG, WGL, NGL, NGSG, SGL, NUM_STATE}; // E=car moving east, W=west, N=sorth, S=south, G=straight|right, L=left turn
	enum {CW_WALK=0, CW_WARN, CW_STOP}; // crosswalk state
	float const state_times[NUM_STATE] = {5.0, 6.0, 5.0, 5.0, 6.0, 5.0}; // in seconds
	unsigned const st_r_orient_masks[NUM_STATE] = {2, 3, 1, 8, 12, 4}; // {W=1, E=2, S=4, N=8}, for straight and right turns
	unsigned const left_orient_masks[NUM_STATE] = {2, 0, 1, 8, 0,  4}; // {W=1, E=2, S=4, N=8}, for left turns only
	unsigned const to_right  [4] = {3, 2, 0, 1}; // {N, S, W, E}
	unsigned const to_left   [4] = {2, 3, 1, 0}; // {S, N, E, W}
	unsigned const other_lane[4] = {1, 0, 3, 2}; // {E, W, N, S}
	unsigned const conn_left[4] = {3,2,0,1}, conn_right[4] = {2,3,1,0};
	colorRGBA const stoplight_colors[3] = {GREEN, YELLOW, RED};
	colorRGBA const crosswalk_colors[3] = {WHITE, ORANGE, ORANGE};

	float stoplight_max_height();

	class stoplight_t {
		uint8_t num_conn, conn, cur_state;
		bool at_conn_road; // longer light times in this case
		float cur_state_ticks;
		// these are mutable because they are set during car update logic, where roads are supposed to be const
		mutable uint8_t car_waiting_sr, car_waiting_left, cw_in_use; // one bit per orient
		mutable bool blocked[4]; // Note: 4 bit flags corresponding to conn bits

		void next_state() {
			++cur_state;
			if (cur_state == NUM_STATE) {cur_state = 0;} // wraparound
		}
		void advance_state();
		bool any_blocked() const {return (blocked[0] || blocked[1] || blocked[2] || blocked[3]);}
		bool is_any_car_waiting_at_this_state() const;
		void find_state_with_waiting_car();
		void run_update_logic();
		float get_cur_state_time_secs() const {return (at_conn_road ? 2.0 : 1.0)*TICKS_PER_SECOND*state_times[cur_state];}
		void ffwd_to_future(float time_secs);
	public:
		stoplight_t(bool at_conn_road_) : num_conn(0), conn(0), cur_state(RED_LIGHT), at_conn_road(at_conn_road_), cur_state_ticks(0.0), car_waiting_sr(0), car_waiting_left(0), cw_in_use(0) {
			reset_blocked();
		}
		void reset_blocked() {UNROLL_4X(blocked[i_] = 0;)}
		void mark_blocked(bool dim, bool dir) const {blocked[2*dim + dir] = 1;} // Note: not actually const, but blocked is mutable
		bool is_blocked(bool dim, bool dir) const {return (blocked[2*dim + dir] != 0);}
		void mark_crosswalk_in_use(bool dim, bool dir) const {cw_in_use |= (1 << (2*dim + dir));}
		void init(uint8_t num_conn_, uint8_t conn_);
		void next_frame();
		void notify_waiting_car(bool dim, bool dir, unsigned turn) const;
		bool red_light(bool dim, bool dir, unsigned turn) const;
		unsigned get_light_state(bool dim, bool dir, unsigned turn) const;
		unsigned get_future_light_state(bool dim, bool dir, unsigned turn, float future_seconds) const;
		bool can_walk(bool dim, bool dir) const;
		unsigned get_crosswalk_state(bool dim, bool dir) const;
		bool check_int_clear(unsigned orient, unsigned turn_dir) const;
		bool check_int_clear(car_base_t const &car) const {return check_int_clear(car.get_orient(), car.turn_dir);}
		bool can_turn_right_on_red(car_base_t const &car) const;
		string str() const;
		string label_str() const;
		colorRGBA get_stoplight_color(bool dim, bool dir, unsigned turn) const {return stoplight_colors[get_light_state(dim, dir, turn)];}
	};
} // end stoplight_ns


struct draw_state_t {
	shader_t s;
	vector3d xlate;
	bool use_building_lights;
	unsigned pass_ix;
protected:
	bool use_smap, use_bmap, shadow_only, use_dlights, emit_now;
	point_sprite_drawer_sized light_psd; // for car/traffic lights
	string label_str;
	point label_pos;
public:
	draw_state_t() : use_building_lights(0), pass_ix(0), use_smap(0), use_bmap(0), shadow_only(0), use_dlights(0), emit_now(0) {}
 	virtual ~draw_state_t() {}
	void set_enable_normal_map(bool val) {use_bmap = val;}
	virtual void draw_unshadowed() {}
	void begin_tile(point const &pos, bool will_emit_now=0, bool ensure_active=0);
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_, bool always_setup_shader);
	void end_draw();
	virtual void post_draw();
	void ensure_shader_active();
	void draw_and_clear_light_flares();
	bool check_sphere_visible(point const &pos, float radius) const {return camera_pdu.sphere_visible_test((pos + xlate), radius);}
	bool check_cube_visible(cube_t const &bc, float dist_scale=1.0, bool shadow_only=0) const;
	static void set_cube_pts(cube_t const &c, float z1f, float z1b, float z2f, float z2b, bool d, bool D, point p[8]);
	static void set_cube_pts(cube_t const &c, float z1, float z2, bool d, bool D, point p[8]) {set_cube_pts(c, z1, z1, z2, z2, d, D, p);}
	static void set_cube_pts(cube_t const &c, bool d, bool D, point p[8]) {set_cube_pts(c, c.z1(), c.z2(), d, D, p);}
	static void rotate_pts(point const &center, float sine_val, float cos_val, int d, int e, point p[8]);
	void draw_cube(quad_batch_draw &qbd, color_wrapper const &cw, point const &center, point const p[8], bool skip_bottom,
		bool invert_normals=0, float tscale=0.0, unsigned skip_dims=0) const;
	void draw_cube(quad_batch_draw &qbd, cube_t const &c, color_wrapper const &cw, bool skip_bottom, float tscale=0.0, unsigned skip_dims=0) const;
	bool add_light_flare(point const &flare_pos, vector3d const &n, colorRGBA const &color, float alpha, float radius);
	void set_label_text(string const &str, point const &pos) {label_str = str; label_pos = pos;}
	void show_label_text();
}; // draw_state_t


namespace streetlight_ns {

	colorRGBA const pole_color(BLACK); // so we don't have to worry about shadows
	colorRGBA const light_color(1.0, 0.9, 0.7, 1.0);
	float const light_height = 0.5; // in units of road width
	float const pole_radius  = 0.015;
	float const light_radius = 0.025;
	float const light_dist   = 3.0;
	float get_streetlight_height();
	float get_streetlight_pole_radius();

	struct streetlight_t {
		point pos; // bottom point
		vector3d dir;
		mutable bool cached_smap;

		streetlight_t(point const &pos_, vector3d const &dir_) : pos(pos_), dir(dir_), cached_smap(0) {}
		bool operator<(streetlight_t const &s) const {return ((pos.y == s.pos.y) ? (pos.x < s.pos.x) : (pos.y < s.pos.y));} // compare y then x
		bool is_lit(bool always_on) const {return (always_on || is_night(STREETLIGHT_ON_RAND*signed_rand_hash(pos.x + pos.y)));}
		point get_lpos() const;
		void draw(draw_state_t &dstate, bool shadow_only, bool is_local_shadow, bool always_on) const;
		void add_dlight(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const;
		bool proc_sphere_coll(point &center, float radius, vector3d const &xlate, vector3d *cnorm) const;
		bool line_intersect(point const &p1, point const &p2, float &t) const;
	};
} // end streetlight_ns


struct streetlights_t {
	vector<streetlight_ns::streetlight_t> streetlights;

	void draw_streetlights(draw_state_t &dstate, bool shadow_only, bool always_on) const;
	void add_streetlight_dlights(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const;
	bool check_streetlight_sphere_coll_xy(point const &center, float radius) const;
	bool proc_streetlight_sphere_coll(point &pos, float radius, vector3d const &xlate, vector3d *cnorm) const;
	bool line_intersect_streetlights(point const &p1, point const &p2, float &t) const;
	void sort_streetlights_by_yx() {sort(streetlights.begin(), streetlights.end());}
};


struct road_isec_t : public cube_t {
	unsigned char num_conn, conn; // connected roads in {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
	short conn_to_city;
	short rix_xy[4], conn_ix[4]; // road/segment index: pos=cur city road, neg=global road; always segment ix
	stoplight_ns::stoplight_t stoplight; // Note: not always needed, maybe should be by pointer/index?

	road_isec_t(cube_t const &c, int rx, int ry, unsigned char conn_, bool at_conn_road, short conn_to_city_=-1);
	tex_range_t get_tex_range(float ar) const;
	void make_4way(unsigned conn_to_city_);
	void next_frame() {stoplight.next_frame();}
	void notify_waiting_car(car_base_t const &car) const {stoplight.notify_waiting_car(car.dim, car.dir, car.turn_dir);}
	bool is_global_conn_int() const {return (rix_xy[0] < 0 || rix_xy[1] < 0 || rix_xy[2] < 0 || rix_xy[3] < 0);}
	bool red_light(car_base_t const &car) const {return stoplight.red_light(car.dim, car.dir, car.turn_dir);}
	bool red_or_yellow_light(car_base_t const &car) const {return (stoplight.get_light_state(car.dim, car.dir, car.turn_dir) != stoplight_ns::GREEN_LIGHT);}
	bool yellow_light(car_base_t const &car) const {return (stoplight.get_light_state(car.dim, car.dir, car.turn_dir) == stoplight_ns::YELLOW_LIGHT);}
	bool will_be_green_light_in(car_base_t const &car, float future_seconds) const {
		return (stoplight.get_future_light_state(car.dim, car.dir, car.turn_dir, future_seconds) == stoplight_ns::GREEN_LIGHT);
	}
	bool can_go_based_on_light(car_base_t const &car) const;
	bool is_orient_currently_valid(unsigned orient, unsigned turn_dir) const;
	unsigned get_dest_orient_for_car_in_isec(car_base_t const &car, bool is_entering) const;
	bool can_go_now(car_t const &car) const;
	bool is_blocked(car_base_t const &car) const {return (can_go_based_on_light(car) && !stoplight.check_int_clear(car));} // light is green but intersection is blocked
	bool has_left_turn_signal(unsigned orient) const;
	cube_t get_stoplight_cube(unsigned n) const;
	bool check_sphere_coll(point const &pos, float radius) const;
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, float dist, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
	void draw_sl_block(quad_batch_draw &qbd, draw_state_t &dstate, point p[4], float h, unsigned state, bool draw_unlit, float flare_alpha, vector3d const &n, tex_range_t const &tr) const;
	void draw_stoplights(quad_batch_draw &qbd, draw_state_t &dstate, bool shadow_only) const;
};


struct road_connector_t : public road_t, public streetlights_t {
	road_t src_road;

	road_connector_t(road_t const &road) : road_t(road), src_road(road) {}
	float get_player_zval(point const &center, cube_t const &c) const;
	void add_streetlights(unsigned num_per_side, bool staggered, float dn_shift_mult, float za, float zb);
};

struct bridge_t : public road_connector_t {
	bool make_bridge;

	bridge_t(road_t const &road) : road_connector_t(road), make_bridge(0) {}
	void add_streetlights() {road_connector_t::add_streetlights(4, 0, 0.05, get_start_z(), get_end_z());} // 4 per side
	bool proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
};

struct tunnel_t : public road_connector_t {
	cube_t ends[2];
	float radius, height, facade_height[2];

	tunnel_t(road_t const &road) : road_connector_t(road), radius(0.0), height(0.0) {facade_height[0] = facade_height[1] = 0.0f;}
	bool enabled() const {return (radius > 0.0);}
	void init(point const &start, point const &end, float radius_, float end_length, bool dim);
	void add_streetlights() {road_connector_t::add_streetlights(2, 1, -0.15, ends[0].z1(), ends[1].z1());} // 2 per side, staggered
	cube_t get_tunnel_bcube() const;
	void calc_top_bot_side_cubes(cube_t cubes[4]) const;
	bool check_mesh_disable(cube_t const &query_region) const {return (ends[0].intersects_xy(query_region) || ends[1].intersects_xy(query_region));} // check both ends
	bool proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const;
	void calc_zvals_and_eext(float &zf, float &zb, float &end_ext) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
};

struct range_pair_t {
	unsigned s, e; // Note: e is one past the end
	range_pair_t(unsigned s_=0, unsigned e_=0) : s(s_), e(e_) {}
	void update(unsigned v);
};

class road_draw_state_t : public draw_state_t {
	quad_batch_draw qbd_batched[NUM_RD_TIDS], qbd_sl, qbd_bridge;
	float ar;

	void draw_city_region_int(quad_batch_draw &cache, unsigned type_ix);
public:
	road_draw_state_t() : ar(1.0) {}
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only, bool always_setup_shader=0);
	virtual void draw_unshadowed();
	virtual void post_draw();
	template<typename T> void add_city_quad(T const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) {add_flat_city_quad(r, qbd, color, ar);} // generic flat road case
	void add_city_quad(road_seg_t  const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) {r.add_road_quad(qbd, color, ar);} // road segment
	void add_city_quad(road_t      const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) {r.add_road_quad(qbd, color, ar/TRACKS_WIDTH);} // tracks
	void add_city_quad(road_plot_t const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool draw_all) { // plots and parks
		if (draw_all || (type_ix == TYPE_PARK) == r.is_park) {add_flat_city_quad(r, qbd, color, ar);}
	}
	template<typename T> void draw_city_region(vector<T> const &v, range_pair_t const &rp, quad_batch_draw &cache, unsigned type_ix, bool draw_all=0) {
		if (rp.s == rp.e) return; // empty
		assert(rp.s <= rp.e);
		assert(rp.e <= v.size());
		assert(type_ix < NUM_RD_TIDS);
		colorRGBA const color(road_colors[type_ix]);

		if (cache.empty()) { // generate and cache quads
			for (unsigned i = rp.s; i < rp.e; ++i) {add_city_quad(v[i], cache, color, type_ix, draw_all);}
		}
		draw_city_region_int(cache, type_ix);
	}
	void draw_bridge(bridge_t const &bridge, bool shadow_only);
	void add_bridge_quad(point const pts[4], color_wrapper const &cw, float normal_scale);
	void draw_tunnel(tunnel_t const &tunnel, bool shadow_only);
	void draw_stoplights(vector<road_isec_t> const &isecs, range_pair_t const &rp, bool shadow_only);
}; // road_draw_state_t

class ao_draw_state_t : public draw_state_t {

protected:
	occlusion_checker_t occlusion_checker;
public:
	quad_batch_draw ao_qbd;
	vect_cube_t &get_occluders() {return occlusion_checker.occluders;}
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_);
	bool is_occluded(cube_t const &bcube) {return (!shadow_only && occlusion_checker.is_occluded(bcube));} // Note: non-const - OC state temp_points is modified
	void draw_ao_qbd();
	virtual void draw_unshadowed() {draw_ao_qbd();}
};

class car_draw_state_t : public ao_draw_state_t { // and trucks and helicopters

	quad_batch_draw qbds[2]; // unshadowed, shadowed
	car_model_loader_t &car_model_loader;
	helicopter_model_loader_t &helicopter_model_loader;
public:
	car_draw_state_t(car_model_loader_t &car_model_loader_, helicopter_model_loader_t &helicopter_model_loader_) :
		car_model_loader(car_model_loader_), helicopter_model_loader(helicopter_model_loader_) {}
	static float get_headlight_dist();
	colorRGBA get_headlight_color(car_t const &car) const;
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_);
	virtual void draw_unshadowed();
	void add_car_headlights(vector<car_t> const &cars, vector3d const &xlate_, cube_t &lights_bcube);
	void gen_car_pts(car_t const &car, bool include_top, point pb[8], point pt[8]) const;
	void draw_car(car_t const &car, bool is_dlight_shadows, bool in_garage);
	void draw_helicopter(helicopter_t const &h, bool shadow_only);
	void add_car_headlights(car_t const &car, cube_t &lights_bcube);
}; // car_draw_state_t


// forward declarations of some classes
class city_road_gen_t;
struct pedestrian_t;
class ped_manager_t;

struct ped_city_vect_t {
	vector<vector<vector<sphere_t>>> peds; // per city per road
	void add_ped(pedestrian_t const &ped, unsigned road_ix);
	void clear();
};

class car_manager_t { // and trucks and helicopters

	car_model_loader_t car_model_loader;
	helicopter_model_loader_t helicopter_model_loader;

	struct car_block_t {
		unsigned start, cur_city, first_parked;
		car_block_t(unsigned s, unsigned c) : start(s), cur_city(c), first_parked(0) {}
		bool is_in_building() const {return (cur_city == NO_CITY_IX);}
	};

	struct helipad_t {
		cube_t bcube;
		bool in_use, reserved;
		helipad_t() : in_use(0), reserved(0) {}
		bool is_avail() const {return !(in_use || reserved);}
	};
	city_road_gen_t const &road_gen;
	vector<car_t> cars;
	vector<car_block_t> car_blocks, car_blocks_by_road;
	vector<cube_with_ix_t> cars_by_road;
	vector<helicopter_t> helicopters;
	vector<helipad_t> helipads;
	ped_city_vect_t peds_crossing_roads;
	car_draw_state_t dstate;
	rand_gen_t rgen;
	vector<unsigned> entering_city;
	cube_t garages_bcube;
	unsigned first_parked_car, first_garage_car;
	bool car_destroyed;

	cube_t get_cb_bcube(car_block_t const &cb ) const;
	road_isec_t const &get_car_isec(car_t const &car) const;
	bool check_collision(car_t &c1, car_t &c2) const;
	void register_car_at_city(car_t const &car);
	cube_t const &get_car_dest_bcube(car_t const &car, bool isec_only) const;
	void add_car();
	void get_car_ix_range_for_cube(vector<car_block_t>::const_iterator cb, cube_t const &bc, unsigned &start, unsigned &end) const;
	void remove_destroyed_cars();
	void update_cars();
	int find_next_car_after_turn(car_t &car);
	void setup_occluders();
	vector3d get_helicopter_size(unsigned model_id);
	void draw_helicopters(bool shadow_only);
public:
	car_manager_t(city_road_gen_t const &road_gen_) :
		road_gen(road_gen_), dstate(car_model_loader, helicopter_model_loader), first_parked_car(0), first_garage_car(0), car_destroyed(0) {}
	bool empty() const {return cars.empty();}
	void clear() {cars.clear(); car_blocks.clear();}
	unsigned get_model_gpu_mem() const {return (car_model_loader.get_gpu_mem() + helicopter_model_loader.get_gpu_mem());}
	void init_cars(unsigned num);
	void add_parked_cars(vector<car_t> const &new_cars, vect_cube_t const &garages);
	void finalize_cars();
	void add_helicopters(vect_cube_t const &hp_locs);
	void extract_car_data(vector<car_city_vect_t> &cars_by_city) const;
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm) const;
	void destroy_cars_in_radius(point const &pos_in, float radius);
	bool get_color_at_xy(point const &pos, colorRGBA &color, int int_ret) const;
	car_t const *get_car_at_pt(point const &pos, bool is_parked) const;
	car_t const *get_car_at(point const &p1, point const &p2) const;
	car_t const *get_car_at_player(float max_dist) const;
	cube_t const &get_car_bcube(unsigned car_id) const {assert(car_id < cars.size()); return cars[car_id].bcube;}
	bool line_intersect_cars(point const &p1, point const &p2, float &t) const;
	bool check_car_for_ped_colls(car_t &car) const;
	void next_frame(ped_manager_t const &ped_manager, float car_speed);
	void helicopters_next_frame(float car_speed);
	bool check_helicopter_coll(cube_t const &bc) const;
	void draw(int trans_op_mask, vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows, bool garages_pass);
	void add_car_headlights(vector3d const &xlate, cube_t &lights_bcube) {dstate.add_car_headlights(cars, xlate, lights_bcube);}
	void free_context() {car_model_loader.free_context(); helicopter_model_loader.free_context();}
}; // car_manager_t


struct pedestrian_t : public waiting_obj_t {

	point target_pos, dest_car_center; // since cars are sorted each frame, we can't find their positions by index so we need to cache them here
	vector3d dir, vel;
	point pos;
	float radius, speed, anim_time, retreat_time;
	unsigned plot, next_plot, dest_plot, dest_bldg; // Note: can probably be made unsigned short later, though these are global plot and building indices
	unsigned short city, model_id, ssn, colliding_ped, cur_rseed;
	unsigned char stuck_count;
	bool collided, ped_coll, is_stopped, in_the_road, at_crosswalk, at_dest, has_dest_bldg, has_dest_car, destroyed, in_building, following_player, is_on_stairs, has_key;

	pedestrian_t(float radius_) : target_pos(all_zeros), dir(zero_vector), vel(zero_vector), pos(all_zeros), radius(radius_), speed(0.0), anim_time(0.0), retreat_time(0.0),
		plot(0), next_plot(0), dest_plot(0), dest_bldg(0), city(0), model_id(0), ssn(0), colliding_ped(0), cur_rseed(1), stuck_count(0), collided(0), ped_coll(0), is_stopped(0),
		in_the_road(0), at_crosswalk(0), at_dest(0), has_dest_bldg(0), has_dest_car(0), destroyed(0), in_building(0), following_player(0), is_on_stairs(0), has_key(0) {}
	bool operator<(pedestrian_t const &ped) const {return ((city == ped.city) ? (plot < ped.plot) : (city < ped.city));} // currently only compares city + plot
	string get_name() const;
	string str() const;
	float get_speed_mult() const;
	float get_height () const {return PED_HEIGHT_SCALE*radius;}
	float get_width  () const {return PED_WIDTH_SCALE *radius;}
	float get_z1     () const {return (pos.z - radius);}
	float get_z2     () const {return (get_z1() + get_height());}
	cube_t get_bcube () const;
	bool target_valid() const {return (target_pos != all_zeros);}
	//bool on_stairs   () const {return (fabs(pos.z - target_pos.z) > 0.01*radius);} // walking on a slope; allow for some floating-point error
	bool on_stairs   () const {return is_on_stairs;}
	bool is_waiting_or_stopped() const {return (speed == 0.0 || waiting_start > 0);}
	void set_velocity(vector3d const &v) {vel = v*(speed/v.mag());} // normalize to original velocity
	void move(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, float &delta_dir);
	void stop();
	void go();
	void wait_for(float seconds);
	bool check_for_safe_road_crossing(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t *dbg_cubes=nullptr) const;
	bool check_ped_ped_coll_range(vector<pedestrian_t> &peds, unsigned pid, unsigned ped_start, unsigned target_plot, float prox_radius, vector3d &force);
	bool check_ped_ped_coll(ped_manager_t const &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, float delta_dir);
	bool check_ped_ped_coll_stopped(vector<pedestrian_t> &peds, unsigned pid);
	bool check_inside_plot(ped_manager_t &ped_mgr, point const &prev_pos, cube_t &plot_bcube, cube_t &next_plot_bcube);
	bool check_road_coll(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube) const;
	bool is_valid_pos(vect_cube_t const &colliders, bool &ped_at_dest, ped_manager_t const *const ped_mgr) const;
	bool try_place_in_plot(cube_t const &plot_cube, vect_cube_t const &colliders, unsigned plot_id, rand_gen_t &rgen);
	point get_dest_pos(cube_t const &plot_bcube, cube_t const &next_plot_bcube, ped_manager_t const &ped_mgr) const;
	bool choose_alt_next_plot(ped_manager_t const &ped_mgr);
	void get_avoid_cubes(ped_manager_t const &ped_mgr, vect_cube_t const &colliders, cube_t const &plot_bcube, cube_t const &next_plot_bcube, point &dest_pos, vect_cube_t &avoid) const;
	void next_frame(ped_manager_t &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, rand_gen_t &rgen, float delta_dir);
	void register_at_dest();
	void destroy() {destroyed = 1;} // that's it, no other effects
	bool is_close_to_player() const;
	void debug_draw(ped_manager_t &ped_mgr) const;
private:
	void run_path_finding(ped_manager_t &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t const &colliders, vector3d &dest_pos);
	void get_plot_bcubes_inc_sidewalks(ped_manager_t const &ped_mgr, cube_t &plot_bcube, cube_t &next_plot_bcube) const;
};

unsigned const MAX_PATH_DEPTH = 32;

class path_finder_t {
	struct path_t : public vector<point> {
		float length;
		path_t() : length(0.0) {}
		path_t(point const &a, point const &b) {init(a, b);} // line constructor
		void init(point const &a, point const &b) {length = p2p_dist(a, b); push_back(a); push_back(b);}
		float calc_length_up_to(const_iterator i) const;
		void calc_length() {length = calc_length_up_to(end());}
		cube_t calc_bcube() const;
	};
	vect_cube_t avoid;
	vector<uint8_t> used;
	path_t path_stack[MAX_PATH_DEPTH];
	float gap;
	point pos, dest;
	cube_t plot_bcube;
	path_t cur_path, best_path, partial_path;
	bool debug;

	bool add_pt_to_path(point const &p, path_t &path) const;
	bool add_pts_around_cube_xy(path_t &path, path_t const &cur_path, path_t::const_iterator p, cube_t const &c, bool dir);
	void find_best_path_recur(path_t const &cur_path, unsigned depth);
	bool shorten_path(path_t &path) const;
public:
	path_finder_t(bool debug_=0) : gap(0.0f), debug(debug_) {}
	vect_cube_t &get_avoid_vector() {return avoid;}
	vector<point> const &get_best_path() const {return (found_complete_path() ? best_path : partial_path);}
	bool found_complete_path() const {return (!best_path.empty());}
	bool found_path() const {return (found_complete_path() || !partial_path.empty());}
	bool find_best_path();
	unsigned run(point const &pos_, point const &dest_, cube_t const &plot_bcube_, float gap_, point &new_dest);
};

class ped_manager_t { // pedestrians

	struct city_ixs_t {
		unsigned ped_ix, plot_ix;
		city_ixs_t() : ped_ix(0), plot_ix(0) {}
		void assign(unsigned ped_ix_, unsigned plot_ix_) {ped_ix = ped_ix_; plot_ix = plot_ix_;}
	};
	city_road_gen_t const &road_gen;
	car_manager_t const &car_manager; // used for ped road crossing safety and dest car selection
	ped_model_loader_t ped_model_loader;
	vector<pedestrian_t> peds, peds_b; // dynamic city, static building
	vector<city_ixs_t> by_city; // first ped/plot index for each city
	vector<unsigned> by_plot;
	vector<unsigned char> need_to_sort_city;
	vector<car_city_vect_t> cars_by_city;
	vector<point> bldg_ppl_pos;
	rand_gen_t rgen;
	ao_draw_state_t dstate;
	int selected_ped_ssn;
	unsigned animation_id;
	bool ped_destroyed, need_to_sort_peds;

	void assign_ped_model(pedestrian_t &ped);
	bool gen_ped_pos(pedestrian_t &ped);
	void expand_cube_for_ped(cube_t &cube) const;
	void remove_destroyed_peds();
	void sort_by_city_and_plot();
	road_isec_t const &get_car_isec(car_base_t const &car) const;
	void register_ped_new_plot(pedestrian_t const &ped);
	int get_road_ix_for_ped_crossing(pedestrian_t const &ped, bool road_dim) const;
	bool draw_ped(pedestrian_t const &ped, shader_t &s, pos_dir_up const &pdu, vector3d const &xlate, float def_draw_dist, float draw_dist_sq,
		bool &in_sphere_draw, bool shadow_only, bool is_dlight_shadows, bool enable_animations);
public:
	// for use in pedestrian_t, mostly for collisions and path finding
	path_finder_t path_finder;
	vect_cube_t const &get_colliders_for_plot(unsigned city_ix, unsigned plot_ix) const;
	road_plot_t const &get_city_plot_for_peds(unsigned city_ix, unsigned plot_ix) const;
	cube_t get_expanded_city_bcube_for_peds(unsigned city_ix) const;
	cube_t get_expanded_city_plot_bcube_for_peds(unsigned city_ix, unsigned plot_ix) const;
	bool is_city_residential(unsigned city_ix) const;
	car_manager_t const &get_car_manager() const {return car_manager;}
	void choose_new_ped_plot_pos(pedestrian_t &ped);
	bool check_isec_sphere_coll(pedestrian_t const &ped) const;
	bool check_streetlight_sphere_coll(pedestrian_t const &ped) const;
	bool mark_crosswalk_in_use(pedestrian_t const &ped);
	bool choose_dest_building_or_parked_car(pedestrian_t &ped);
	unsigned get_next_plot(pedestrian_t &ped, int exclude_plot=-1) const;
	void move_ped_to_next_plot(pedestrian_t &ped);
	bool has_nearby_car(pedestrian_t const &ped, bool road_dim, float delta_time, vect_cube_t *dbg_cubes=nullptr) const;
	bool has_nearby_car_on_road(pedestrian_t const &ped, bool dim, unsigned road_ix, float delta_time, vect_cube_t *dbg_cubes) const;
	bool has_car_at_pt(point const &pos, unsigned city, bool is_parked) const;
	bool choose_dest_parked_car(unsigned city_id, unsigned &plot_id, unsigned &car_ix, point &car_center);
public:
	ped_manager_t(city_road_gen_t const &road_gen_, car_manager_t const &car_manager_) :
		road_gen(road_gen_), car_manager(car_manager_), selected_ped_ssn(-1), animation_id(1), ped_destroyed(0), need_to_sort_peds(0) {}
	void next_animation();
	static float get_ped_radius();
	bool empty() const {return (peds.empty() && peds_b.empty());}
	void clear() {peds.clear(); peds_b.clear(); by_city.clear();}
	unsigned get_model_gpu_mem() const {return ped_model_loader.get_gpu_mem();}
	void init(unsigned num_city, unsigned num_building);
	bool proc_sphere_coll(point &pos, float radius, vector3d *cnorm) const;
	bool line_intersect_peds(point const &p1, point const &p2, float &t) const;
	void destroy_peds_in_radius(point const &pos_in, float radius);
	void next_frame();
	pedestrian_t const *get_ped_at(point const &p1, point const &p2) const;
	unsigned get_first_ped_at_plot(unsigned plot) const {assert(plot < by_plot.size()); return by_plot[plot];}
	void get_peds_crossing_roads(ped_city_vect_t &pcv) const;
	void draw(vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows);
	void draw_peds_in_building(int first_ped_ix, ped_draw_vars_t const &pdv);
	void get_ped_bcubes_for_building(int first_ped_ix, unsigned bix, vect_cube_t &bcubes, bool moving_only) const;
	void register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity);
	void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only);
	void free_context() {ped_model_loader.free_context();}
	//vector3d get_dest_move_dir(point const &pos) const;
}; // end ped_manager_t


template <typename T> void remove_destroyed(vector<T> &objs) {
	typename vector<T>::iterator i(objs.begin()), o(i);
	for (; i != objs.end(); ++i) {if (!i->destroyed) {*(o++) = *i;}}
	objs.erase(o, objs.end());
}

bool check_line_clip_update_t(point const &p1, point const &p2, float &t, cube_t const &c);
point rand_xy_pt_in_cube(cube_t const &c, float radius, rand_gen_t &rgen);
bool sphere_in_light_cone_approx(pos_dir_up const &pdu, point const &center, float radius);

// from gen_buildings.cpp
bool enable_building_people_ai();
bool place_building_people(vect_building_place_t &locs, float radius, float speed_mult, unsigned num);
void update_building_ai_state(vector<pedestrian_t> &people, float delta_dir);
void get_all_garages(vect_cube_t &garages);
void get_all_city_helipads(vect_cube_t &helipads);
bool check_city_building_line_coll_bs(point const &p1, point const &p2, point &p_int);
void update_buildings_zmax_for_line(point const &p1, point const &p2, float radius, float house_extra_zval, float &cur_zmax);
