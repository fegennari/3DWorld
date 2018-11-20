// 3D World - City Header
// by Frank Gennari
// 11/19/18

#pragma once

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"
#include "model3d.h"

using std::string;

unsigned  const CONN_CITY_IX((1<<16)-1); // uint16_t max

enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   TID_PARK_LOT,  TID_TRACKS,  NUM_RD_TIDS};
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, TYPE_PARK_LOT, TYPE_TRACKS, NUM_RD_TYPES};
enum {TURN_NONE=0, TURN_LEFT, TURN_RIGHT, TURN_UNSPEC};
enum {INT_NONE=0, INT_ROAD, INT_PLOT, INT_PARKING};
enum {RTYPE_ROAD=0, RTYPE_TRACKS};
unsigned const CONN_TYPE_NONE = 0;
colorRGBA const road_colors[NUM_RD_TYPES] = {WHITE, WHITE, WHITE, WHITE, WHITE, WHITE, WHITE}; // parking lots are darker than roads

int       const FORCE_MODEL_ID = -1; // -1 disables
unsigned  const NUM_CAR_COLORS = 10;
colorRGBA const car_colors[NUM_CAR_COLORS] = {WHITE, GRAY_BLACK, GRAY, ORANGE, RED, DK_RED, DK_BLUE, DK_GREEN, YELLOW, BROWN};

float const CONN_ROAD_SPEED_MULT = 2.0; // twice the speed limit on connector roads
float const HEADLIGHT_ON_RAND    = 0.1;
vector3d const CAR_SIZE(0.30, 0.13, 0.08); // {length, width, height} in units of road width

extern double tfticks;
extern float fticks;


inline bool is_isect(unsigned type) {return (type >= TYPE_ISEC2 && type <= TYPE_ISEC4);}
inline int encode_neg_ix(unsigned ix) {return -(int(ix)+1);}
inline unsigned decode_neg_ix(int ix) {assert(ix < 0); return -(ix+1);}

inline float rand_hash(float to_hash) {return fract(12345.6789*to_hash);}
inline float signed_rand_hash(float to_hash) {return 0.5*(rand_hash(to_hash) - 1.0);}

struct car_model_t {

	string fn;
	int body_mat_id, fixed_color_id;
	float xy_rot, dz, lod_mult, scale; // xy_rot in degrees
	vector<unsigned> shadow_mat_ids;

	car_model_t() : body_mat_id(-1), fixed_color_id(-1), xy_rot(0.0), dz(0.0), lod_mult(1.0), scale(1.0) {}
	car_model_t(string const &fn_, int bmid, int fcid, float rot, float dz_, float lm, vector<unsigned> const &smids) :
		fn(fn_), body_mat_id(bmid), fixed_color_id(fcid), xy_rot(rot), dz(dz_), lod_mult(lm), shadow_mat_ids(smids) {}
	bool read(FILE *fp);
};

struct city_params_t {

	unsigned num_cities, num_samples, num_conn_tries, city_size_min, city_size_max, city_border, road_border, slope_width, num_rr_tracks;
	float road_width, road_spacing, conn_road_seg_len, max_road_slope;
	unsigned make_4_way_ints; // 0=all 3-way intersections; 1=allow 4-way; 2=all connector roads must have at least a 4-way on one end; 4=only 4-way (no straight roads)
							  // cars
	unsigned num_cars;
	float car_speed, traffic_balance_val, new_city_prob, max_car_scale;
	bool enable_car_path_finding;
	vector<car_model_t> car_model_files;
	// parking lots
	unsigned min_park_spaces, min_park_rows;
	float min_park_density, max_park_density;
	// lighting
	bool car_shadows;
	unsigned max_lights, max_shadow_maps;
	// trees
	unsigned max_trees_per_plot;
	float tree_spacing;
	// detail objects
	unsigned max_benches_per_plot;

	city_params_t() : num_cities(0), num_samples(100), num_conn_tries(50), city_size_min(0), city_size_max(0), city_border(0), road_border(0), slope_width(0),
		num_rr_tracks(0), road_width(0.0), road_spacing(0.0), conn_road_seg_len(1000.0), max_road_slope(1.0), make_4_way_ints(0), num_cars(0), car_speed(0.0),
		traffic_balance_val(0.5), new_city_prob(1.0), max_car_scale(1.0), enable_car_path_finding(0), min_park_spaces(12), min_park_rows(1), min_park_density(0.0),
		max_park_density(1.0), car_shadows(0), max_lights(1024), max_shadow_maps(0), max_trees_per_plot(0), tree_spacing(1.0), max_benches_per_plot(0) {}
	bool enabled() const {return (num_cities > 0 && city_size_min > 0);}
	bool roads_enabled() const {return (road_width > 0.0 && road_spacing > 0.0);}
	float get_road_ar() const {return nearbyint(road_spacing/road_width);} // round to nearest texture multiple
	static bool read_error(string const &str) {cout << "Error reading city config option " << str << "." << endl; return 0;}
	bool read_option(FILE *fp);
	vector3d get_nom_car_size() const {return CAR_SIZE*road_width;}
	vector3d get_max_car_size() const {return max_car_scale*get_nom_car_size();}
}; // city_params_t


struct car_t;

struct road_gen_base_t {
	virtual cube_t get_bcube_for_car(car_t const &car) const = 0;
};


struct car_t {
	cube_t bcube, prev_bcube;
	bool dim, dir, stopped_at_light, entering_city, in_tunnel, dest_valid, destroyed;
	unsigned char cur_road_type, color_id, turn_dir, front_car_turn_dir, model_id;
	unsigned short cur_city, cur_road, cur_seg, dest_city, dest_isec;
	float height, dz, rot_z, turn_val, cur_speed, max_speed, waiting_pos, waiting_start;
	car_t const *car_in_front;

	car_t() : bcube(all_zeros), dim(0), dir(0), stopped_at_light(0), entering_city(0), in_tunnel(0), dest_valid(0), destroyed(0), cur_road_type(TYPE_RSEG),
		color_id(0), turn_dir(TURN_NONE), front_car_turn_dir(TURN_UNSPEC), model_id(0), cur_city(0), cur_road(0), cur_seg(0), dest_city(0), dest_isec(0),
		height(0.0), dz(0.0), rot_z(0.0), turn_val(0.0), cur_speed(0.0), max_speed(0.0), waiting_pos(0.0), waiting_start(0.0), car_in_front(nullptr) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	point get_center() const {return bcube.get_cube_center();}
	unsigned get_orient() const {return (2*dim + dir);}
	unsigned get_orient_in_isec() const {return (2*dim + (!dir));} // invert dir (incoming, not outgoing)
	float get_max_speed() const {return ((cur_city == CONN_CITY_IX) ? CONN_ROAD_SPEED_MULT : 1.0)*max_speed;}
	float get_length() const {return (bcube.d[ dim][1] - bcube.d[ dim][0]);}
	float get_width () const {return (bcube.d[!dim][1] - bcube.d[!dim][0]);}
	float get_max_lookahead_dist() const;
	bool is_almost_stopped() const {return (cur_speed < 0.1*max_speed);}
	bool is_stopped () const {return (cur_speed == 0.0);}
	bool is_parked  () const {return (max_speed == 0.0);}
	bool in_isect   () const {return is_isect(cur_road_type);}
	bool headlights_on() const {return (!is_parked() && (in_tunnel || is_night(HEADLIGHT_ON_RAND*signed_rand_hash(height + max_speed))));} // no headlights when parked
	unsigned get_isec_type() const {assert(in_isect()); return (cur_road_type - TYPE_ISEC2);}
	void park() {cur_speed = max_speed = 0.0;}
	float get_turn_rot_z(float dist_to_turn) const;
	float get_wait_time_secs  () const {return (float(tfticks) - waiting_start)/TICKS_PER_SECOND;} // Note: only meaningful for cars stopped at lights
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
	void stop() {cur_speed = 0.0;} // immediate stop
	void move_by(float val) {bcube.d[dim][0] += val; bcube.d[dim][1] += val;}
	bool check_collision(car_t &c, road_gen_base_t const &road_gen);
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d const &xlate, vector3d *cnorm) const;
	point get_front(float dval=0.5) const;
	bool front_intersects_car(car_t const &c) const;
	void honk_horn_if_close() const;
	void honk_horn_if_close_and_fast() const;
	void on_alternate_turn_dir(rand_gen_t &rgen);
	void register_adj_car(car_t &c);
	unsigned count_cars_in_front(cube_t const &range=cube_t(all_zeros)) const;
	float get_sum_len_space_for_cars_in_front(cube_t const &range) const;
};


struct comp_car_road_then_pos {
	vector3d const &xlate;
	comp_car_road_then_pos(vector3d const &xlate_) : xlate(xlate_) {}
	bool operator()(car_t const &c1, car_t const &c2) const;
};


class car_model_loader_t : public model3ds {
	vector<int> models_valid;
	void ensure_models_loaded() {if (empty()) {load_car_models();}}
public:
	static unsigned num_models();

	bool is_model_valid(unsigned id);
	car_model_t const &get_model(unsigned id) const;
	void load_car_models();
	void draw_car(shader_t &s, vector3d const &pos, cube_t const &car_bcube, vector3d const &dir, colorRGBA const &color,
		point const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail);
};

