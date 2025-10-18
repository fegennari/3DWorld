// 3D World - City Header
// by Frank Gennari
// 11/19/18

#pragma once

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"
#include "shaders.h"
#include "draw_utils.h"
#include "buildings.h"
#include "city_model.h"

using std::string;

unsigned const CONN_CITY_IX((1<<16)-1); // uint16_t max
unsigned const NO_CITY_IX(CONN_CITY_IX-1); // used for cars not in any city (in house garages)
unsigned const CITY_BIX_START((1<<15)); // start of building rooftop indices

// Note: when addint to this list, must also update road_mat_mgr_t::ensure_road_textures()
enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   TID_PARK_LOT,  TID_TRACKS,  TID_PARK,  TID_DRIVEWAY,  TID_ROAD_SKIRT,  TID_TURN_SKIRT, /*none for bldg*/ NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, TYPE_PARK_LOT, TYPE_TRACKS, TYPE_PARK, TYPE_DRIVEWAY, TYPE_ROAD_SKIRT, TYPE_TURN_SKIRT, TYPE_BUILDING,   NUM_RD_TYPES};
enum {TURN_NONE=0, TURN_LEFT, TURN_RIGHT, TURN_UNSPEC};
enum {INT_NONE=0, INT_ROAD, INT_PLOT, INT_PARKING, INT_PARK, INT_TRACK, INT_BUILDING, INT_TURBINE};
enum {RTYPE_ROAD=0, RTYPE_TRACKS};
unsigned const CONN_TYPE_NONE = 0;
colorRGBA const road_colors[NUM_RD_TYPES] = {WHITE, WHITE, WHITE, WHITE, WHITE, WHITE, WHITE, LT_GRAY, LT_GRAY, GRAY, GRAY, WHITE}; // white except for parks, driveways, and skirts

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
float const SIGN_STOPSIGN_HEIGHT = 2.2; // height of street sign relative to attached stop sign
vector3d const CAR_SIZE(0.30, 0.13, 0.08); // {length, width, height} in units of road width
float const CAR_RADIUS_SCALE(CAR_SIZE.mag()/CAR_SIZE.z);


inline bool is_isect(unsigned type) {return (type >= TYPE_ISEC2 && type <= TYPE_ISEC4);}
inline int encode_neg_ix(unsigned ix) {return -(int(ix)+1);}
inline unsigned decode_neg_ix(int ix) {assert(ix < 0); return -(ix+1);}

inline float rand_hash(float to_hash) {return fract(12345.6789*to_hash);}
inline float signed_rand_hash(float to_hash) {return 0.5*(rand_hash(to_hash) - 1.0);}


struct city_params_t {

	unsigned num_cities=0, num_samples=100, num_conn_tries=50, city_size_min=0, city_size_max=0, city_border=0, road_border=0, slope_width=0, num_rr_tracks=0, park_rate=0;
	float road_width=0.0, road_spacing=0.0, road_spacing_rand=0.0, road_spacing_xy_add=0.0, conn_road_seg_len=1000.0, max_road_slope=1.0, max_track_slope=1.0;
	float residential_probability=0.0, model_anim_scale=1.0;
	unsigned make_4_way_ints=0; // 0=all 3-way intersections; 1=allow 4-way; 2=all connector roads must have at least a 4-way on one end; 3=only 4-way (no straight roads)
	unsigned add_tlines=2; // 0=never, 1=always, 2=only when there are no secondary buildings
	bool assign_house_plots=0, new_city_conn_road_alg=0, add_skyways=0;
	// cars
	unsigned num_cars=0;
	float car_speed=0.0, traffic_balance_val=0.5, new_city_prob=1.0, max_car_scale=1.0;
	bool enable_car_path_finding=0, convert_model_files=0, cars_use_driveways=0;
	vector<city_model_t> car_model_files, ped_model_files, hc_model_files;
	// parking lots
	unsigned min_park_spaces=12, min_park_rows=1;
	float min_park_density=0.0, max_park_density=1.0;
	// lighting
	bool car_shadows=0;
	unsigned max_lights=1024, max_shadow_maps=0, smap_size=0;
	// trees
	unsigned max_trees_per_plot=0;
	float tree_spacing=1.0;
	// detail objects
	unsigned max_benches_per_plot=0;
	// pedestrians
	unsigned num_peds=0;
	float ped_speed=0.0;
	bool ped_respawn_at_dest=0, use_animated_people=0;
	bool any_model_has_animations=0; // calculated, not specified in the config file
	string default_anim_name;
	// buildings; maybe should be building params, but we have the model loading code here
	vector<city_model_t> building_models[NUM_OBJ_MODELS]; // multiple model files per type
	// use for option reading
	int read_error_flag=0;
	kw_to_val_map_t<bool     >  kwmb;
	kw_to_val_map_t<unsigned >  kwmu;
	kw_to_val_map_float_check_t kwmr;

	city_params_t() : kwmb(read_error_flag, "city"), kwmu(read_error_flag, "city"), kwmr(read_error_flag, "city") {init_kw_maps();}
	bool enabled() const {return (num_cities > 0 && city_size_min > 0);}
	bool roads_enabled() const {return (road_width > 0.0 && road_spacing > 0.0);}
	float get_road_ar () const {return round(road_spacing/road_width);} // round to nearest texture multiple
	static bool read_error(string const &str) {cout << "Error reading city config option " << str << "." << endl; return 0;}
	bool read_option(FILE *fp);
	bool add_model(unsigned id, FILE *fp);
	bool has_helicopter_model() const {return !hc_model_files.empty();}
	vector3d get_nom_car_size() const {return CAR_SIZE*road_width;}
	vector3d get_max_car_size() const {return max_car_scale*get_nom_car_size();}
private:
	void init_kw_maps();
}; // city_params_t


struct car_base_t;
struct driveway_t;
struct wind_turbine_t;
class road_draw_state_t;

struct road_gen_base_t {
	virtual cube_t get_bcube_for_car(car_base_t const &car) const = 0;
 	virtual ~road_gen_base_t() {}
};

struct car_base_t { // the part needed for the pedestrian interface (size = 48); could use oriented_cube_t
	cube_t bcube;
	bool dim=0, dir=0, stopped_at_light=0, stopped_for_ssign=0, need_gas=0; // Note: stopped_at_light also applies to stopped at stop sign
	uint8_t cur_road_type=TYPE_RSEG, turn_dir=TURN_NONE;
	uint16_t cur_city=0, cur_road=0, cur_seg=0;
	short dest_driveway=-1, dest_gstation=-1; // -1 is unset
	float max_speed=0.0, cur_speed=0.0;

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
	float get_front_end() const {return bcube.d[dim][ dir];}
	float get_back_end () const {return bcube.d[dim][!dir];}
	point get_front(float dval=0.5) const;
};

struct car_t : public car_base_t, public waiting_obj_t { // size = 136
	bool is_truck=0, is_police=0, is_ambulance=0, is_emergency=0, entering_city=0, in_tunnel=0, dest_valid=0, destroyed=0;
	bool in_reverse=0, engine_running=0, is_braking=0, in_parking_lot=0;
	uint8_t color_id=0, front_car_turn_dir=TURN_UNSPEC, model_id=0;
	uint16_t dest_city=0, dest_isec=0, dest_gs_lane=0;
	float height=0.0, dz=0.0, rot_z=0.0, turn_val=0.0, waiting_pos=0.0, wake_time=0.0;
	vector2d park_space_cent; // or gas station pos
	cube_t prev_bcube;
	car_t const *car_in_front=nullptr;

	void set_bcube(point const &center, vector3d const &sz);
	bool is_valid   () const {return !bcube.is_all_zeros();}
	bool is_sleeping() const {return (wake_time > 0.0);}
	float get_max_lookahead_dist() const;
	bool headlights_on() const;
	float get_turn_rot_z(float dist_to_turn) const;
	bool is_close_to_player() const;
	colorRGBA const &get_color() const;
	unsigned get_unique_id() const {return (unsigned(1000000.0*max_speed) + color_id + (model_id<<8));} // not guaranteed to be unique, but pretty close
	void apply_scale(float scale);
	void set_correct_len_width_from_model(vector3d const &model_sz);
	void destroy();
	float get_min_sep_dist_to_car(car_t const &c, bool add_one_car_len=0) const;
	string str() const;
	string label_str() const;
	void choose_max_speed(rand_gen_t &rgen);
	void move(float speed_mult);
	void set_target_speed(float speed_factor);
	void maybe_accelerate(float mult=0.02);
	void accelerate(float mult=0.02);
	void decelerate(float mult=0.05);
	void decelerate_fast();
	void park() {cur_speed = max_speed = 0.0;}
	void stop() {cur_speed = 0.0;} // immediate stop
	void sleep(rand_gen_t &rgen, float min_time_secs);
	bool maybe_wake(rand_gen_t &rgen);
	void move_by(float val) {bcube.translate_dim(dim, val);}
	void begin_turn() {turn_val = bcube.get_center_dim(!dim);}
	bool maybe_apply_turn(float centerline, bool for_driveway);
	void complete_turn_and_swap_dim();
	void person_in_the_way(bool is_player, bool at_stopsign);
	bool must_wait_entering_or_crossing_road(vector<car_t> const &cars, driveway_t const &driveway, unsigned road_ix, float lookahead_time) const;
	bool check_for_road_clear_and_wait(vector<car_t> const &cars, driveway_t const &driveway, unsigned road_ix);
	// driveway/parking lot/gas station enter/exit
	bool run_enter_driveway_logic(vector<car_t> const &cars, driveway_t const &driveway);
	void pull_into_driveway(driveway_t const &driveway, rand_gen_t &rgen);
	void back_or_pull_out_of_driveway(driveway_t const &driveway);
	bool exit_driveway_to_road(vector<car_t> const &cars, driveway_t const &driveway, float centerline, unsigned road_ix, rand_gen_t &rgen);
	// collisions
	bool check_collision(car_t &c, road_gen_base_t const &road_gen);
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm) const;
	bool front_intersects_car(car_t const &c) const;
	void honk_horn_if_close() const;
	void honk_horn_if_close_and_fast() const;
	void on_alternate_turn_dir(rand_gen_t &rgen);
	void register_adj_car(car_t &c);
	unsigned count_cars_in_front(cube_t const &range=cube_t(all_zeros)) const;
	float get_sum_len_space_for_cars_in_front(cube_t const &range) const;
	cube_t get_parking_space_debug_marker() const;
	cube_t get_ped_coll_check_area() const;
};

struct car_city_vect_t {
	vector<car_base_t> cars[2][2]; // {dim x dir}
	vect_cube_with_ix_t parked_car_bcubes, sleeping_car_bcubes, parking_lot_car_bcubes; // stores car bcube + road/plot/driveway index
	void clear_cars();
};


struct comp_car_road {
	bool operator()(car_base_t const &c1, car_base_t const &c2) const {return (c1.cur_road < c2.cur_road);}
};
struct comp_car_city_then_road {
	bool operator()(car_base_t const &c1, car_base_t const &c2) const {
		if (c1.cur_city != c2.cur_city) return (c1.cur_city < c2.cur_city);
		return ((c1.is_parked() != c2.is_parked()) ? c2.is_parked() : (c1.cur_road < c2.cur_road));
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
	float wait_time=0.0; // time to wait before takeoff
	float fly_zval=0.0; // zval required for flight to avoid buildings and terrain
	float blade_rot=0.0;
	unsigned dest_hp=0; // destination (or current) helipad
	unsigned state=STATE_WAIT, model_id=0;
	bool dynamic=0, dynamic_shadow=0;

	helicopter_t(cube_t const &bcube_, vector3d const &dir_, unsigned model_id_, unsigned dest_hp_, bool dynamic_) :
		bcube(bcube_), dir(dir_), dest_hp(dest_hp_), model_id(model_id_), dynamic(dynamic_) {}
	point get_landing_pt() const {return cube_bot_center(bcube);}
	void invalidate_tile_shadow_map(vector3d const &shadow_offset, bool repeat_next_frame) const;
};


class road_mat_mgr_t {

	bool inited=0;
	unsigned tids[NUM_RD_TIDS] = {}, sl_tid=0; // stoplight tid
public:
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
	unsigned road_ix=0;
	float bt_lo=0.0, bt_hi=0.0; // start and end pos of bridge or tunnel in dim
	//unsigned char type; // road, railroad, etc. {RTYPE_ROAD, RTYPE_TRACKS}
	bool dim=0; // dim the road runs in
	bool slope=0; // 0: z1 applies to first (lower) point; 1: z1 applies to second (upper) point
	bool has_bridge=0, has_tunnel=0;

	road_t(cube_t const &c, bool dim_, bool slope_=0, unsigned road_ix_=0) : cube_t(c), road_ix(road_ix_), dim(dim_), slope(slope_) {}
	road_t(point const &s, point const &e, float width, bool dim_, bool slope_=0, unsigned road_ix_=0);
	road_t() {} // only used for name generation for city connector roads
	float get_length   () const {return get_sz_dim( dim);}
	float get_width    () const {return get_sz_dim(!dim);}
	float get_slope_val() const {return dz()/get_length();}
	float get_start_z  () const {return (slope ? z2() : z1());}
	float get_end_z    () const {return (slope ? z1() : z2());}
	float get_z_adj    () const {return (ROAD_HEIGHT + 0.5*get_slope_val()*(dim ? DY_VAL : DX_VAL));} // account for a half texel of error for sloped roads
	void register_bridge_or_tunnel(cube_t const &bc, bool is_bridge);
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, -ar, (dim ? -1.0 : 1.0), 0, dim);}
	cube_t const &get_bcube() const {return *this;}
	cube_t       &get_bcube()       {return *this;}
	void add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar, bool add_skirt=0) const;
	string get_name(unsigned city_ix) const;
};

struct road_seg_t : public road_t {
	unsigned short road_ix, conn_ix[2], conn_type[2]; // {dim=0, dim=1}
	mutable unsigned short car_count=0; // can be written to during car update logic

	void init_ixs() {conn_ix[0] = conn_ix[1] = 0; conn_type[0] = conn_type[1] = CONN_TYPE_NONE;}
	road_seg_t(road_t const &r, unsigned rix) : road_t(r), road_ix(rix) {init_ixs();}
	road_seg_t(cube_t const &c, unsigned rix, bool dim_, bool slope_=0) : road_t(c, dim_, slope_), road_ix(rix) {init_ixs();}
	void next_frame() {car_count = 0;}
};

struct driveway_t : public oriented_cube_t {
	// dim/dir: direction to road; d[dim][dir] is the edge shared with the road
	// in_use is modified by car_manager in a different thread - must be mutable, maybe should be atomic
	mutable uint8_t in_use=0; // either reserves the spot, or a car is parked there; 1=temporary, 2=permanent
	mutable unsigned last_ped_frame=0;
	unsigned plot_ix=0;
	int park_lot_ix=-1, gstation_ix=-1; // driveway may be part of a parking lot or gas station
	float stop_loc=0.0; // used for gas stations

	driveway_t() {}
	driveway_t(cube_t const &c, bool dim_, bool dir_, unsigned pix, int plix=-1, int gsix=-1, float sl=0.0) :
		oriented_cube_t(c, dim_, dir_), plot_ix(pix), park_lot_ix(plix), gstation_ix(gsix), stop_loc(sl) {}
	float get_edge_at_road() const {return d[dim][dir];}
	float get_centerline  () const {return get_center_dim(!dim);}
	void mark_ped_this_frame() const;
	bool has_recent_ped() const;
	bool is_parking_lot() const {return (park_lot_ix >= 0);}
	bool is_gas_station() const {return (gstation_ix >= 0);}
	cube_t extend_across_road() const;
	tex_range_t get_tex_range(float ar) const;
};

struct dw_query_t {
	driveway_t const *const driveway;
	unsigned dix;
	dw_query_t() : driveway(nullptr), dix(0) {}
	dw_query_t(driveway_t const *const driveway_, unsigned dix_) : driveway(driveway_), dix(dix_) {}
};

struct road_plot_t : public cube_t {
	uint8_t xpos, ypos; // position within the city grid
	bool is_residential=0, has_parking=0, is_park=0;
	road_plot_t(cube_t const &c, uint8_t xpos_, uint8_t ypos_, bool is_res=0) : cube_t(c), xpos(xpos_), ypos(ypos_), is_residential(is_res) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, ar, ar);}
	bool is_residential_not_park() const {return ( is_residential && !is_park);}
	bool is_commercial          () const {return (!is_residential && !is_park);}
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

struct parking_lot_t : public oriented_cube_t {
	bool row_dir;
	unsigned short row_sz, num_rows, orig_ix;
	vector<unsigned char> used_spaces;
	parking_lot_t(cube_t const &c, bool dim_, bool dir_, bool rdir, unsigned row_sz_, unsigned num_rows_, unsigned ix) :
		oriented_cube_t(c, dim_, dir_), row_dir(rdir), row_sz(row_sz_), num_rows(num_rows_), orig_ix(ix) {}
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
		uint8_t num_conn=0, conn=0, cur_state=RED_LIGHT;
		bool at_conn_road; // longer light times in this case
		float cur_state_ticks=0.0;
		// these are mutable because they are set during car update logic, where roads are supposed to be const
		mutable uint8_t car_waiting_sr=0, car_waiting_left=0, cw_in_use=0; // one bit per orient
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
		stoplight_t(bool at_conn_road_) : at_conn_road(at_conn_road_) {reset_blocked();}
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

class hedge_draw_t : public vao_manager_t {
	unsigned num_verts=0;
	cube_t bcube;
	vect_cube_t to_draw;

	void create(cube_t const &bc);
public:
	bool empty() const {return to_draw.empty();}
	void add(cube_t const &bc) {to_draw.push_back(bc);}
	void draw_and_clear(shader_t &s);
};

struct draw_state_t {
	shader_t s;
	vector3d xlate;
	point camera_bs;
	bool use_building_lights=0;
	unsigned pass_ix=0;
	float draw_tile_dist=0.0;
	hedge_draw_t hedge_draw;
	vector<vert_wrap_t> temp_verts; // used for sphere drawing
protected:
	bool use_smap=0, use_bmap=0, shadow_only=0, use_dlights=0, emit_now=0;
	point_sprite_drawer_sized light_psd; // for car/traffic lights
	occlusion_checker_t occlusion_checker;
	string label_str;
	point label_pos;
public:
 	virtual ~draw_state_t() {}
	void copy_from(draw_state_t &d) {xlate = d.xlate; camera_bs = d.camera_bs; draw_tile_dist = d.draw_tile_dist; use_smap = d.use_smap; shadow_only = d.shadow_only;}
	void set_enable_normal_map(bool val) {use_bmap = val;}
	bool normal_maps_enabled() const {return use_bmap;}
	vect_cube_t &get_occluders() {return occlusion_checker.occluders;}
	bool is_visible_and_unoccluded(cube_t const &c, float dist_scale=1.0) const {return (check_cube_visible(c, dist_scale) && !is_occluded(c));}
	bool is_occluded(cube_t const &c) const {return (!shadow_only && occlusion_checker.is_occluded(c));}
	float get_lod_factor(point const &pos) const {return draw_tile_dist/p2p_dist(camera_bs, pos);}
	virtual bool has_unshadowed () const {return 0;}
	virtual void draw_unshadowed() {}
	void begin_tile(point const &pos, bool will_emit_now=0, bool ensure_active=0);
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_, bool always_setup_shader,
		bool enable_animations=0, bool enable_occlusion=1, bool enable_reflect=0);
	void end_draw();
	virtual void post_draw();
	void set_untextured_material();
	void unset_untextured_material();
	void ensure_shader_active();
	void draw_and_clear_light_flares();
	bool check_sphere_visible(point const &pos, float radius) const {return camera_pdu.sphere_visible_test((pos + xlate), radius);}
	bool check_cube_visible(cube_t const &bc, float dist_scale=1.0) const;
	static void set_cube_pts(cube_t const &c, float z1f, float z1b, float z2f, float z2b, bool d, bool D, point p[8]);
	static void set_cube_pts(cube_t const &c, float z1, float z2, bool d, bool D, point p[8]) {set_cube_pts(c, z1, z1, z2, z2, d, D, p);}
	static void set_cube_pts(cube_t const &c, bool d, bool D, point p[8]) {set_cube_pts(c, c.z1(), c.z2(), d, D, p);}
	static void rotate_pts(point const &center, float sine_val, float cos_val, int d, int e, point p[8]);
	void draw_cube(quad_batch_draw &qbd, color_wrapper const &cw, point const &center, point const p[8], bool skip_bottom,
		bool invert_normals=0, float tscale=0.0, unsigned skip_dims=0) const;
	void draw_cube(quad_batch_draw &qbd, cube_t const &c, color_wrapper const &cw, bool skip_bottom=0, float tscale=0.0, unsigned skip_dims=0,
		bool mirror_x=0, bool mirror_y=0, bool swap_tc_xy=0, float tscale_x=1.0, float tscale_y=1.0, float tscale_z=1.0, bool skip_top=0, bool no_cull=0) const;
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
		unsigned plot_ix;
		bool on_bridge_or_tunnel=0;
		mutable bool cached_smap=0;

		streetlight_t(point const &pos_, vector3d const &dir_, bool bt=0, unsigned plot_ix_=0) : pos(pos_), dir(dir_), plot_ix(plot_ix_), on_bridge_or_tunnel(bt) {}
		bool operator<(streetlight_t const &s) const {return ((pos.y == s.pos.y) ? (pos.x < s.pos.x) : (pos.y < s.pos.y));} // compare y then x
		bool is_lit(bool always_on) const {return (always_on || is_night(STREETLIGHT_ON_RAND*signed_rand_hash(pos.x + pos.y)));}
		point get_lpos() const;
		void draw(road_draw_state_t &dstate, bool shadow_only, bool is_local_shadow, bool always_on) const;
		void add_dlight(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const;
		bool proc_sphere_coll(point &center, float radius, vector3d const &xlate, vector3d *cnorm) const;
		bool line_intersect(point const &p1, point const &p2, float &t) const;
	};
} // end streetlight_ns


struct streetlights_t {
	vector<streetlight_ns::streetlight_t> streetlights;

	void draw_streetlights(road_draw_state_t &dstate, bool shadow_only, bool always_on) const;
	void add_streetlight_dlights(vector3d const &xlate, cube_t &lights_bcube, bool always_on) const;
	bool check_streetlight_sphere_coll_xy(point const &center, float radius, cube_t &coll_cube) const;
	bool proc_streetlight_sphere_coll(point &pos, float radius, vector3d const &xlate, vector3d *cnorm) const;
	bool line_intersect_streetlights(point const &p1, point const &p2, float &t) const;
	void sort_streetlights_by_yx() {sort(streetlights.begin(), streetlights.end());}
};


struct ssign_state_t { // per incoming orient
	bool in_use=0, is_truck=0, is_emergency=0;
	uint8_t turn_dir=0, dest_orient=0;
	int arrive_frame=0;
};
struct ssign_state_pair_t {
	ssign_state_t waiting, entering;
};

struct road_isec_t : public cube_t {
	bool has_stoplight=0, has_stopsign=0;
	unsigned char num_conn, conn; // connected roads in {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
	unsigned char hospital_dir=0; // 8 bit flags in dim pairs for each side
	short conn_to_city;
	short rix_xy[4], conn_ix[4]; // road/segment index: pos=cur city road, neg=global road; always segment ix
	stoplight_ns::stoplight_t stoplight; // Note: not always needed, maybe should be by pointer/index?
	mutable ssign_state_pair_t ssign_state[4]; // one per orient; updated by cars through const functions, must be mutable

	road_isec_t(cube_t const &c, int rx, int ry, unsigned char conn_, bool at_conn_road, bool has_stoplight_, short conn_to_city_=-1);
	tex_range_t get_tex_range(float ar) const;
	void init_stoplights();
	void make_stop_signs();
	void make_4way(unsigned conn_to_city_);
	void next_frame();
	void notify_waiting_car(car_t const &car) const;
	void notify_leaving_car(car_t const &car) const;
	void notify_turned_car (car_t const &car) const;
	void mark_crosswalk_in_use(bool dim, bool dir) const;
	bool is_global_conn_int() const {return (rix_xy[0] < 0 || rix_xy[1] < 0 || rix_xy[2] < 0 || rix_xy[3] < 0);}
	bool red_light(car_base_t const &car) const {return (has_stoplight && stoplight.red_light(car.dim, car.dir, car.turn_dir));}
	bool red_or_yellow_light(car_base_t const &car) const {return (has_stoplight && stoplight.get_light_state(car.dim, car.dir, car.turn_dir) != stoplight_ns::GREEN_LIGHT);}
	bool yellow_light(car_base_t const &car) const {return (has_stoplight && stoplight.get_light_state(car.dim, car.dir, car.turn_dir) == stoplight_ns::YELLOW_LIGHT);}
	bool will_be_green_light_in(car_base_t const &car, float future_seconds) const {
		return (has_stoplight && stoplight.get_future_light_state(car.dim, car.dir, car.turn_dir, future_seconds) == stoplight_ns::GREEN_LIGHT);
	}
	bool can_go_based_on_light(car_base_t const &car) const;
	bool is_orient_currently_valid(unsigned orient, unsigned turn_dir) const;
	unsigned get_dest_orient_for_car_in_isec(car_base_t const &car, bool is_entering=1) const;
	bool can_go_now(car_t const &car) const;
	bool is_blocked(car_base_t const &car) const {return (can_go_based_on_light(car) && !stoplight.check_int_clear(car));} // light is green but intersection is blocked
	bool has_left_turn_signal(unsigned orient) const;
	cube_t get_stoplight_cube(unsigned n) const;
	point get_stop_sign_pos  (unsigned n) const;
	float get_stop_sign_height  () const {return 0.075*(dx() + dy());}
	float get_street_sign_height() const {return SIGN_STOPSIGN_HEIGHT*get_stop_sign_height();}
	bool check_sphere_coll(point const &pos, float radius, cube_t &coll_cube) const;
	void add_stoplight_bcubes_in_region(cube_t const &region, vect_cube_t &bcubes) const;
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, float dist, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
	void draw_stoplights_and_street_signs(road_draw_state_t &dstate, vector<road_t> const &roads, unsigned cur_city, bool shadow_only) const;
	void add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar, bool add_skirt=0) const;
private:
	void draw_sl_block(quad_batch_draw &qbd, draw_state_t &dstate, point p[4], float h, unsigned state,
		bool draw_unlit, float flare_alpha, vector3d const &n, tex_range_t const &tr) const;
	ssign_state_pair_t &get_ssign_state(car_t const &car) const {return ssign_state[car.get_orient_in_isec()];}
	void init_ssign_state(car_t const &car, ssign_state_t &ss, bool is_entering) const;
};


struct road_connector_t : public road_t, public streetlights_t {
	road_t src_road;

	road_connector_t(road_t const &road) : road_t(road), src_road(road) {}
	float get_player_zval(point const &center, cube_t const &c) const;
	void add_streetlights(unsigned num_per_side, bool staggered, float dn_shift_mult, float za, float zb);
};

struct bridge_t : public road_connector_t {
	bool make_bridge=0, over_water=0;
	float zmin_below=0.0; // for terrain below the bridge

	bridge_t(road_t const &road) : road_connector_t(road) {}
	cube_t get_drawn_bcube() const;
	void add_streetlights() {road_connector_t::add_streetlights(4, 0, 0.05, get_start_z(), get_end_z());} // 4 per side
	bool proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
};

struct tunnel_t : public road_connector_t {
	cube_t ends[2];
	float radius=0.0, height=0.0, facade_height[2]={};

	tunnel_t(road_t const &road) : road_connector_t(road) {}
	bool enabled() const {return (radius > 0.0);}
	void init(point const &start, point const &end, float radius_, float end_length, bool dim);
	void add_streetlights() {road_connector_t::add_streetlights(2, 1, -0.15, ends[0].z1(), ends[1].z1());} // 2 per side, staggered
	cube_t get_tunnel_bcube() const;
	void calc_top_bot_side_cubes(cube_t cubes[4]) const;
	bool check_mesh_disable(cube_t const &query_region) const {return (ends[0].intersects_xy(query_region) || ends[1].intersects_xy(query_region));} // check both ends
	bool proc_sphere_coll(point &center, point const &prev, float sradius, float prev_frame_zval, vector3d const &xlate, vector3d *cnorm) const;
	float get_end_ext() const;
	void calc_zvals_and_eext(float &zf, float &zb, float &end_ext) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
};

struct transmission_line_t {
	unsigned city1, city2;
	float tower_height;
	point p1, p2, p1_wire_pts[3], p2_wire_pts[3];
	cube_t bcube;
	vector<point> tower_pts;
	vector<pair<point, point>> connections;

	transmission_line_t(unsigned c1, unsigned c2, float tower_height_, point const &p1_, point const &p2_) :
		city1(c1), city2(c2), tower_height(tower_height_), p1(p1_), p2(p2_) {}
	void calc_bcube();
	bool sphere_intersect_xy(point const &pos, float radius) const;
	bool cube_intersect_xy(cube_t const &c) const;
};

struct range_pair_t {
	unsigned s, e; // Note: e is one past the end
	range_pair_t(unsigned s_=0, unsigned e_=0) : s(s_), e(e_) {}
	void clear() {s = e = 0;}
	void update(unsigned v);
};

class road_draw_state_t : public draw_state_t {
	quad_batch_draw qbd_batched[NUM_RD_TIDS], qbd_bridge;
public: // used directly by stoplight drawing
	// {stoplight, untextured, streetlight emissive spot, road skirts, hospital sign}; could add qbd_ssign here for stop signs
	quad_batch_draw qbd_sl, qbd_untextured, qbd_emissive, qbd_skirt, qbd_hospital;
	vector<vert_norm_comp_tc_color> text_verts;
	vect_cube_t plot_cuts;
private:
	float ar=1.0;

	void draw_city_region_int(quad_batch_draw &cache, unsigned type_ix);
	void draw_transmission_line_wires(point const &p1, point const &p2, point const wire_pts1[3], point const wire_pts2[3], float radius);
public:
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only, bool always_setup_shader=0, bool enable_occlusion=1, bool enable_reflect=0);
	virtual bool has_unshadowed () const;
	virtual void draw_unshadowed();
	virtual void post_draw();
	void end_cur_tile();
	void add_city_quad(road_seg_t  const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool);
	void add_city_quad(road_isec_t const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool);
	void add_city_quad(road_t      const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool);
	void add_city_quad(road_plot_t const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool draw_all);
	void draw_city_skirt(cube_t const &bcube, bool shadow_only);

	template<typename T> void add_city_quad(T const &r, quad_batch_draw &qbd, colorRGBA const &color, unsigned type_ix, bool) { // generic flat road case
		add_flat_city_quad(r, qbd, color, ar);
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
	void draw_stoplights_and_street_signs(vector<road_isec_t> const &isecs, vector<road_t> const &roads, range_pair_t const &rp, unsigned cur_city, bool shadow_only);
	void draw_transmission_line(transmission_line_t const &tline);
}; // road_draw_state_t

void draw_and_clear_blur_qbd(quad_batch_draw &qbd);

class ao_draw_state_t : public draw_state_t {
public:
	quad_batch_draw ao_qbd;
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_, bool enable_animations=0);
	virtual bool has_unshadowed () const {return !ao_qbd.empty();}
	virtual void draw_unshadowed() {draw_and_clear_blur_qbd(ao_qbd);}
};

class car_draw_state_t : public ao_draw_state_t { // and trucks and helicopters

	quad_batch_draw qbds[2]; // unshadowed, shadowed
	car_model_loader_t &car_model_loader;
	helicopter_model_loader_t &helicopter_model_loader;
	uint64_t last_smap_tile_id=0;
	bool last_smap_tile_id_valid=0;
public:
	car_draw_state_t(car_model_loader_t &car_model_loader_, helicopter_model_loader_t &helicopter_model_loader_) :
		car_model_loader(car_model_loader_), helicopter_model_loader(helicopter_model_loader_) {}
	static float get_headlight_dist();
	colorRGBA get_headlight_color(car_t const &car) const;
	void pre_draw(vector3d const &xlate_, bool use_dlights_, bool shadow_only_);
	virtual bool has_unshadowed () const {return (!qbds[0].empty() || ao_draw_state_t::has_unshadowed());}
	virtual void draw_unshadowed();
	void draw_remaining_cars();
	void add_car_headlights(vector<car_t> const &cars, vector3d const &xlate_, cube_t &lights_bcube);
	static void gen_car_pts(car_t const &car, bool include_top, point pb[8], point pt[8]);
	void draw_car(car_t const &car, bool is_dlight_shadows);
	void draw_helicopter(helicopter_t const &h, bool shadow_only);
	void add_car_headlights(car_t const &car, cube_t &lights_bcube) const;
}; // car_draw_state_t


// forward declarations of some classes
class city_road_gen_t;
class ped_manager_t;
class city_spectate_manager_t;

struct pedestrian_t : public person_base_t { // city pedestrian
	point dest_car_center; // since cars are sorted each frame, we can't find their positions by index so we need to cache them here
	unsigned plot=0, next_plot=0, dest_plot=0, dest_bldg=0; // Note: can probably be made unsigned short later, though these are global plot and building indices
	unsigned short city=0, colliding_ped=0;
	unsigned char stuck_count=0, reverse_count=0;
	bool collided=0, ped_coll=0, in_the_road=0, at_crosswalk=0, at_dest=0, has_dest_bldg=0, has_dest_car=0, destroyed=0, follow_player=0, using_nav_grid=0;

	pedestrian_t(float radius_) : person_base_t(radius_) {}
	bool operator<(pedestrian_t const &ped) const {return ((city == ped.city) ? (plot < ped.plot) : (city < ped.city));} // currently only compares city + plot
	std::string str() const;
	float get_coll_radius() const {return 0.6f*radius;} // using a smaller radius to allow peds to get close to each other
	float get_speed_mult () const;
	void destroy() {destroyed = 1;} // that's it, no other effects
	void clear_current_dest() {has_dest_bldg = has_dest_car = at_dest = 0;}
	void move(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, float &delta_dir);
	void update_velocity_dir(vector3d const &force, float delta_dir);
	bool check_for_safe_road_crossing(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t *dbg_cubes=nullptr) const;
	bool check_ped_ped_coll_range(vector<pedestrian_t> &peds, unsigned pid, unsigned ped_start, unsigned target_plot, float prox_radius, vector3d &force);
	void run_collision_avoid(point const &ipos, vector3d const &ivel, float r2, float dist_sq, bool is_player, vector3d &force);
	bool overlaps_player_in_z(point const &player_pos) const;
	bool check_ped_ped_coll(ped_manager_t const &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, float delta_dir);
	bool check_ped_ped_coll_stopped(vector<pedestrian_t> &peds, unsigned pid);
	bool check_inside_plot(ped_manager_t &ped_mgr, point const &prev_pos, cube_t &plot_bcube, cube_t &next_plot_bcube);
	bool check_road_coll(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, cube_t &coll_cube) const;
	bool is_valid_pos(vect_cube_t const &colliders, bool &ped_at_dest, cube_t &coll_cube, ped_manager_t const *const ped_mgr, int *coll_bldg_ix=nullptr) const;
	bool is_possibly_valid_dest_pos(point const &dpos, vect_cube_t const &colliders) const;
	bool try_place_in_plot(cube_t const &plot_cube, vect_cube_t const &colliders, unsigned plot_id, rand_gen_t &rgen);
	point get_dest_pos(cube_t const &plot_bcube, cube_t const &next_plot_bcube, ped_manager_t const &ped_mgr, int &debug_state) const;
	bool choose_alt_next_plot(ped_manager_t const &ped_mgr);
	void get_avoid_cubes(ped_manager_t const &ped_mgr, vect_cube_t const &colliders, cube_t const &plot_bcube, cube_t const &next_plot_bcube,
		point &dest_pos, vect_cube_t &avoid, bool &in_illegal_area, bool &avoid_entire_plot) const;
	bool check_path_blocked(ped_manager_t &ped_mgr, point const &dest, bool check_buildings);
	void next_frame(ped_manager_t &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, rand_gen_t &rgen, float delta_dir);
	void register_at_dest();
	void debug_draw(ped_manager_t &ped_mgr) const;
private:
	bool can_target_player(ped_manager_t const &ped_mgr) const;
	void run_path_finding(ped_manager_t &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t const &colliders, vector3d &dest_pos);
	void get_plot_bcubes_inc_sidewalks(ped_manager_t const &ped_mgr, cube_t &plot_bcube, cube_t &next_plot_bcube) const;
};

struct peds_per_city_t {
	vector<vect_sphere_t> by_road, by_p_lot;
};
struct ped_city_vect_t {
	vector<peds_per_city_t> peds; // per city
	void add_ped(pedestrian_t const &ped, unsigned rp_ix, bool in_parking_lot);
	void clear();
};

class car_manager_t { // and trucks and helicopters

	car_model_loader_t car_model_loader;
	helicopter_model_loader_t helicopter_model_loader;

	struct car_block_t {
		unsigned start, cur_city, first_parked=0;
		cube_t bcube;
		car_block_t(unsigned s, unsigned c, cube_t const &bc) : start(s), cur_city(c), bcube(bc) {}
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
	vector<cube_with_ix_t> cars_by_road; // Note: car_blocks_by_road and cars_by_road are used by map mode
	vector<helicopter_t> helicopters;
	vector<helipad_t> helipads;
	ped_city_vect_t peds_crossing_roads;
	car_draw_state_t dstate;
	rand_gen_t rgen;
	vector<unsigned> entering_city;
	unsigned first_parked_car=0;
	bool car_destroyed=0;

	road_isec_t const &get_car_isec(car_t const &car) const;
	bool check_collision(car_t &c1, car_t &c2) const;
	void register_car_at_city(car_t const &car);
	cube_t get_car_dest_bcube(car_t const &car, bool isec_only) const;
	int get_parking_lot_ix_for_car(car_t const &car) const;
	void add_car();
	void get_car_ix_range_for_cube(vector<car_block_t>::const_iterator cb, cube_t const &bc, unsigned &start, unsigned &end) const;
	void remove_destroyed_cars();
	void update_cars();
	int find_next_car_after_turn(car_t &car);
	void setup_occluders();
	vector3d get_helicopter_size(unsigned model_id);
	vector<bridge_t      > const &get_bridges      () const;
	vector<wind_turbine_t> const &get_wind_turbines() const;
	void draw_helicopters(bool shadow_only);
public:
	friend class city_spectate_manager_t;
	car_manager_t(city_road_gen_t const &road_gen_) : road_gen(road_gen_), dstate(car_model_loader, helicopter_model_loader) {}
	bool empty() const {return cars.empty();}
	bool has_car_models() const {return !car_model_loader.empty();}
	size_t get_model_gpu_mem() const {return (car_model_loader.get_gpu_mem() + helicopter_model_loader.get_gpu_mem());}
	void clear();
	void init_cars(unsigned num);
	void add_parked_cars(vector<car_t> const &new_cars);
	void finalize_cars();
	void assign_car_model_size_color(car_t &car, rand_gen_t &local_rgen, bool is_in_garage, unsigned btype=BTYPE_UNSET);
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
	void draw(int trans_op_mask, vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows);
	void set_car_model_color(car_t &car, unsigned btype);
	void draw_car_in_pspace(car_t &car, shader_t &s, vector3d const &xlate, bool shadow_only, unsigned btype);
	void add_car_headlights(vector3d const &xlate, cube_t &lights_bcube) {dstate.add_car_headlights(cars, xlate, lights_bcube);}
	void free_context() {car_model_loader.free_context(); helicopter_model_loader.free_context();}
}; // car_manager_t


unsigned const MAX_PATH_DEPTH = 32;

class path_finder_t {
	struct path_t : public vector<point> {
		float length=0.0;
		path_t() {}
		path_t(point const &a, point const &b) {init(a, b);} // line constructor
		void init(point const &a, point const &b) {length = p2p_dist(a, b); push_back(a); push_back(b);}
		float calc_length_up_to(const_iterator i) const;
		void calc_length() {length = calc_length_up_to(end());}
		cube_t calc_bcube() const;
	};
	vect_cube_t avoid;
	vector<uint8_t> used;
	path_t path_stack[MAX_PATH_DEPTH];
	float gap=0.0;
	point pos, dest, prev_target_pos;
	cube_t plot_bcube;
	path_t cur_path, best_path, partial_path;
	bool debug=0;

	bool add_pt_to_path(point const &p, path_t &path) const;
	bool add_pts_around_cube_xy(path_t &path, path_t const &cur_path, path_t::const_iterator p, cube_t const &c, bool dir);
	void find_best_path_recur(path_t const &cur_path, unsigned depth);
	bool shorten_path(path_t &path) const;
public:
	path_finder_t(bool debug_=0) : debug(debug_) {}
	vect_cube_t &get_avoid_vector() {return avoid;}
	vector<point> const &get_best_path() const {return (found_complete_path() ? best_path : partial_path);}
	bool found_complete_path() const {return (!best_path.empty());}
	bool found_path() const {return (found_complete_path() || !partial_path.empty());}
	bool find_best_path();
	unsigned run(point const &pos_, point const &dest_, point const &prev_target_pos_, cube_t const &plot_bcube_, float gap_, point &new_dest);
};

class city_cube_nav_grid_manager;

class ped_manager_t { // pedestrians

	struct city_ixs_t {
		unsigned ped_ix, plot_ix;
		city_ixs_t() : ped_ix(0), plot_ix(0) {}
		void assign(unsigned ped_ix_, unsigned plot_ix_) {ped_ix = ped_ix_; plot_ix = plot_ix_;}
	};
	city_road_gen_t const &road_gen;
	car_manager_t const &car_manager; // used for ped road crossing safety and dest car selection
	ped_model_loader_t ped_model_loader;
	vector<pedestrian_t> peds; // dynamic city pedestrians
	vector<city_ixs_t> by_city; // first ped/plot index for each city
	vector<unsigned> by_plot;
	vector<unsigned char> need_to_sort_city;
	car_city_vect_t empty_cars_vect;
	vector<car_city_vect_t> cars_by_city;
	vector<person_t const *> to_draw;
	rand_gen_t rgen;
	ao_draw_state_t dstate;
	unique_ptr<city_cube_nav_grid_manager> nav_grid_mgr;
	int selected_ped_ssn=-1;
	unsigned animation_id=ANIM_ID_WALK, tot_num_plots=0;
	bool ped_destroyed=0, need_to_sort_peds=0, prev_choose_zombie=0;

	void assign_ped_model(person_base_t &ped);
	void maybe_reassign_ped_model(person_base_t &ped);
	bool gen_ped_pos(pedestrian_t &ped);
	void expand_cube_for_ped(cube_t &cube) const;
	void remove_destroyed_peds();
	void sort_by_city_and_plot();
	road_isec_t const &get_car_isec(car_base_t const &car) const;
	int get_road_ix_for_ped_crossing(pedestrian_t const &ped, bool road_dim     ) const;
	int get_parking_lot_ix_for_ped  (pedestrian_t const &ped, bool inc_driveways) const;
	void setup_occluders();
	bool draw_ped(person_base_t const &ped, shader_t &s, pos_dir_up const &pdu, vector3d const &xlate, float def_draw_dist, float draw_dist_sq,
		bool &in_sphere_draw, bool shadow_only, bool is_dlight_shadows, animation_state_t *anim_state, bool is_in_building);
	car_city_vect_t const &get_cars_for_city(unsigned city) const {return ((city < cars_by_city.size()) ? cars_by_city[city] : empty_cars_vect);}
public:
	friend class city_spectate_manager_t;
	// for use in pedestrian_t, mostly for collisions and path finding
	path_finder_t path_finder;
	ai_path_t grid_path;

	ped_manager_t(city_road_gen_t const &road_gen_, car_manager_t const &car_manager_);
	ped_manager_t (ped_manager_t const &) = delete; // forbidden
	void operator=(ped_manager_t const &) = delete; // forbidden
	~ped_manager_t();
	city_cube_nav_grid_manager &get_nav_grid_mgr();
	vect_cube_t const &get_colliders_for_plot(unsigned city_ix, unsigned plot_ix) const;
	road_plot_t const &get_city_plot_for_peds(unsigned city_ix, unsigned plot_ix) const;
	int get_global_plot_id_for_pos(unsigned city_ix, point const &pos) const;
	void register_ped_new_plot(pedestrian_t const &ped);
	dw_query_t get_nearby_driveway(unsigned city_ix, unsigned global_plot_ix, point const &pos, float dist) const;
	car_base_t const *find_car_using_driveway(unsigned city_ix, dw_query_t const &dw) const;
	cube_t get_city_bcube_for_peds(unsigned city_ix) const;
	cube_t get_expanded_city_bcube_for_peds(unsigned city_ix) const;
	cube_t get_expanded_city_plot_bcube_for_peds(unsigned city_ix, unsigned plot_ix) const;
	bool is_city_residential(unsigned city_ix) const;
	car_manager_t const &get_car_manager() const {return car_manager;}
	void choose_new_ped_plot_pos(pedestrian_t &ped);
	bool check_isec_sphere_coll       (pedestrian_t const &ped, cube_t &coll_cube) const;
	bool check_streetlight_sphere_coll(pedestrian_t const &ped, cube_t &coll_cube) const;
	bool mark_crosswalk_in_use(pedestrian_t const &ped);
	bool choose_dest_building_or_parked_car(pedestrian_t &ped);
	unsigned get_tot_num_plots() const {return tot_num_plots;}
	unsigned get_next_plot(pedestrian_t &ped, int exclude_plot=-1) const;
	void move_ped_to_next_plot(pedestrian_t &ped);
	// cars
	bool has_cars_in_city(unsigned city_ix) const {return (city_ix < cars_by_city.size());}
	bool has_nearby_car(pedestrian_t const &ped, bool road_dim, float delta_time, vect_cube_t *dbg_cubes=nullptr) const;
	bool has_nearby_car_on_road(pedestrian_t const &ped, bool dim, unsigned road_ix, float delta_time, vect_cube_t *dbg_cubes) const;
	bool has_car_at_pt(point const &pos, unsigned city, bool is_parked) const;
	bool has_parked_car_on_path(point const &p1, point const &p2, unsigned city) const;
	void get_parked_car_bcubes_for_plot(cube_t const &plot, unsigned city, vect_cube_t &car_bcubes) const;
	bool choose_dest_parked_car(unsigned city_id, unsigned &plot_id, unsigned &car_ix, point &car_center);
	void next_animation();
	static float get_ped_radius();
	void clear();
	void init(unsigned num_city);
	size_t get_model_gpu_mem() const {return ped_model_loader.get_gpu_mem();}
	void maybe_reassign_models();
	bool proc_sphere_coll(point &pos, float radius, vector3d *cnorm) const;
	bool line_intersect_peds(point const &p1, point const &p2, float &t) const;
	void destroy_peds_in_radius(point const &pos_in, float radius);
	void next_frame();
	pedestrian_t const *get_ped_at(point const &p1, point const &p2) const;
	unsigned get_first_ped_at_plot(unsigned plot) const {assert(plot < by_plot.size()); return by_plot[plot];}
	void get_peds_crossing_roads(ped_city_vect_t &pcv) const;
	void get_pedestrians_in_area(cube_t const &area, int building_ix, vector<point> &pts) const;
	void draw(vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows);
	void gen_and_draw_people_in_building(ped_draw_vars_t const &pdv);
	void draw_people_in_building(vector<person_t> const &people, ped_draw_vars_t const &pdv);
	unsigned get_player_model_id();
	bool is_player_model_female();
	void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only);
	void free_context() {ped_model_loader.free_context();}
}; // end ped_manager_t


template <typename T> void remove_destroyed(vector<T> &objs) {
	typename vector<T>::iterator i(objs.begin()), o(i);
	for (; i != objs.end(); ++i) {if (!i->destroyed) {*(o++) = *i;}}
	objs.erase(o, objs.end());
}

bool check_line_clip_update_t(point const &p1, point const &p2, float &t, cube_t const &c);
point rand_xy_pt_in_cube(cube_t const &c, float radius, rand_gen_t &rgen);
bool sphere_in_light_cone_approx(pos_dir_up const &pdu, point const &center, float radius);
bool moving_sphere_cube_intersect_xy(point const &p1, point const &p2, cube_t const &c, float dist, float radius);
cube_t get_city_bcube(unsigned city_id);
cube_t get_city_bcube_at_pt(point const &pos);
float get_sidewalk_width();
float get_inner_sidewalk_width();
void        set_z_plane_rect_pts  (point const &center, float rx, float ry, point pts[4]);
inline void set_z_plane_square_pts(point const &center, float radius,       point pts[4]) {set_z_plane_rect_pts(center, radius, radius, pts);}
point get_player_pos_bs();
// from gen_buildings.cpp
bool have_city_buildings();
bool enable_building_people_ai();
void update_building_ai_state(float delta_dir);
void get_all_city_helipads(vect_cube_t &helipads);
cube_t get_cur_basement();
bool check_city_building_line_coll_bs(point const &p1, point const &p2, point &p_int);
void update_buildings_zmax_for_line(point const &p1, point const &p2, float radius, float house_extra_zval, float &cur_zmax);
bool check_sphere_coll_building(point const &pos, float radius, bool xy_only, unsigned building_id);
// from city_interact.cpp
void init_city_spectate_manager(car_manager_t &car_manager, ped_manager_t &ped_manager);
bool skip_bai_draw(person_t     const &bai);
bool skip_ped_draw(pedestrian_t const &ped);
bool skip_car_draw(car_t        const &car);

