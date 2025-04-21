// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h" // for vbo_wrap_t
#include "transform_obj.h" // for xform_matrix
#include "draw_utils.h" // for quad_batch_draw
#include "file_utils.h" // for kw_to_val_map_t
#include "pedestrians.h"
#include "building_animals.h"

bool const EXACT_MULT_FLOOR_HEIGHT = 1;
bool const ENABLE_MIRROR_REFLECTIONS = 1;
bool const ENABLE_GLASS_FLOOR_REF  = 1;
bool const DRAW_CITY_INT_WINDOWS   = 1; // requires having different window x/y size/space/offset values for interior vs. exterior windows
bool const ADD_WALKWAY_EXT_DOORS   = 1; // requires DRAW_WALKWAY_INTERIORS=1 in gen_buildings.cpp
unsigned const MAX_CYLIN_SIDES     = 36;
unsigned const MAX_DRAW_BLOCKS     = 8; // for building interiors only; currently have floor, ceiling, walls, and doors
unsigned const NUM_STAIRS_PER_FLOOR= 12;
unsigned const NUM_STAIRS_PER_FLOOR_U = 16;
unsigned const NUM_STAIRS_PER_FLOOR_L = 12;
unsigned const NUM_STAIRS_PER_FLOOR_ESC = 12;
unsigned const MAX_SHELVES         = 4;
float const FLOOR_THICK_VAL_HOUSE  = 0.10; // 10% of floor spacing
float const FLOOR_THICK_VAL_OFFICE = 0.11; // thicker for office buildings
float const FLOOR_THICK_VAL_WINDOWLESS = 0.12; // even thicker for windowless office buildings
float const RAMP_THICKNESS_SCALE   = 0.11; // thickness to height ratio
float const WINDOW_BORDER_MULT     = 0.94; // account for the frame part of the window texture, which is included in the interior cutout of the window
float const WALL_THICK_VAL         = 0.05; // 5% of floor spacing
float const DOOR_THICK_TO_WIDTH    = 0.04; // ratio of door thickness to width for doors opening to the side
float const DEF_CITY_MIN_ALPHA     = 0.01;
float const DOOR_WIDTH_SCALE       = 0.5;
float const DOOR_WIDTH_SCALE_OFFICE= 0.7; // wider than house doors
float const STAIRS_WALL_WIDTH_MULT = 0.15; // relative to the depth of a stair
float const ELEVATOR_Z2_SHIFT      = 0.6; // shift downward, relative to ceiling thickness
float const ELEVATOR_STAND_DIST    = 0.4; // distance that people stand from the elevator door relative to floor spacing
float const ELEVATOR_FRAME_WIDTH   = 0.2; // frame width to each side of door relative to elevator width
float const DOOR_FRAME_WIDTH       = 0.07; // for door texture, relative to door width
float const EXT_BASEMENT_JOIN_DIST = 4.0; // relative to floor spacing
float const BALCONY_PILLAR_SCALE   = 0.15; // relative to depth
float const BASEMENT_ENTRANCE_SCALE= 0.33;
float const CHAIR_LEG_WIDTH        = 0.15; // relative to chair width
float const CHAIR_LEG_WIDTH_MALL   = 0.07; // relative to chair width
float const BED_HEAD_WIDTH         = 0.04; // headboard width; relative to bed width
float const SHELF_RACK_HEIGHT_FS   = 0.85*(1.0 - FLOOR_THICK_VAL_OFFICE);
float const WHEATER_PIPE_SPACING   = 0.65; // relative to radius
float const WHEATER_PIPE_H_DIST    = 0.92; // relative to radius
float const PLANT_POT_RADIUS       = 0.75; // relative to plant bcube radius
float const DEF_NORM_BIAS_SCALE    = 10.0; // see shadow_map.part
float const CITY_BIAS_SCALE        = 0.1;
float const RETAIL_SMAP_DSCALE     = 0.7; // player shadow caster distance for retail rooms and malls relative to light radius
float const PERSON_INT_SMAP_DSCALE = 0.8; // building person shadow caster distance relative to light radius
float const ESCALATOR_SPEED        = 1.0/TICKS_PER_SECOND; // in steps per tick
float const MALL_FLOOR_HEIGHT      = 2.0; // as a multiple of normal building floor height
float const CEILING_BEAM_THICK     = 2.5; // as a multiple of wall thickness
float const BACKSPLASH_HEIGHT      = 0.33; // relative to cabinet height

unsigned const NUM_CHAIR_COLORS = 13;
unsigned const MAX_BCASE_BOOKS  = 48; // limited by available bit flags
unsigned const NUM_BOOK_COLORS  = 16;
unsigned const NUM_PAPER_COLORS = 6;
unsigned const NUM_SPCAN_COLORS = 11;
unsigned const NUM_LAMP_COLORS  = 6;
unsigned const NUM_TCAN_COLORS  = 6;
unsigned const NUM_TAPE_COLORS  = 7;
unsigned const NUM_SHIRT_COLORS = 14;
unsigned const NUM_STAPLER_COLORS = 5;
unsigned const NUM_TSHIRT_COLORS  = 9;
unsigned const NUM_TOASTER_COLORS = 7;
unsigned const NUM_TRASH_COLORS   = 8;
unsigned const NUM_PFLOAT_COLORS  = 6;
unsigned const NUM_SP_EMISSIVE_COLORS = 2;
colorRGBA const GD_SP_COLOR(0.5, 1.0, 1.0); // used for glow-in-the-dark spraypaint
colorRGBA const chair_colors[NUM_CHAIR_COLORS] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN, DK_BROWN};
colorRGBA const book_colors [NUM_BOOK_COLORS ] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, DK_BLUE, BLUE, LT_BLUE, DK_RED, RED, ORANGE, YELLOW, DK_GREEN, LT_BROWN, BROWN, DK_BROWN};
colorRGBA const spcan_colors[NUM_SPCAN_COLORS] = {GD_SP_COLOR, WHITE, RED, GREEN, BLUE, YELLOW, PINK, ORANGE, PURPLE, BROWN, BLACK};
colorRGBA const sp_emissive_colors[NUM_SP_EMISSIVE_COLORS] = {colorRGBA(0.2, 1.0, 0.2), colorRGBA(0.2, 0.5, 1.0)}; // light green, greenish blue
colorRGBA const lamp_colors[NUM_LAMP_COLORS]   = {WHITE, GRAY_BLACK, BROWN, LT_BROWN, DK_BROWN, OLIVE};
colorRGBA const cream(0.9, 0.9, 0.8), vlt_yellow(1.0, 1.0, 0.5);
colorRGBA const paper_colors[NUM_PAPER_COLORS] = {WHITE, WHITE, WHITE, cream, cream, vlt_yellow};
colorRGBA const pen_colors    [4] = {WHITE, BLACK, colorRGBA(0.2, 0.4, 1.0), RED};
colorRGBA const pencil_colors [2] = {colorRGBA(1.0, 0.75, 0.25), colorRGBA(1.0, 0.5, 0.1)};
colorRGBA const marker_colors [8] = {BLACK, RED, BLACK, BLUE, BLACK, GREEN, RED, PURPLE};
colorRGBA const railing_colors[3] = {GOLD, DK_BROWN, BLACK};
colorRGBA const tcan_colors   [NUM_TCAN_COLORS   ] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
colorRGBA const tape_colors   [NUM_TAPE_COLORS   ] = {GRAY, GRAY, GRAY, GRAY, BKGRAY, BKGRAY, colorRGBA(0.2, 0.2, 1.0)}; // gray duct tape is the most common
colorRGBA const shirt_colors  [NUM_SHIRT_COLORS  ] = {WHITE, WHITE, WHITE, BKGRAY, GRAY, RED, DK_RED, LT_BLUE, BLUE, DK_BLUE, DK_GREEN, DK_BROWN, BROWN, ORANGE};
colorRGBA const stapler_colors[NUM_STAPLER_COLORS] = {BLACK, RED, BLACK, BLUE, BLACK};
colorRGBA const TSHIRT_COLORS [NUM_TSHIRT_COLORS ] = {WHITE, BKGRAY, GRAY, RED, GREEN, BLUE, YELLOW, ORANGE, WHITE};
colorRGBA const toaster_colors[NUM_TOASTER_COLORS] = {WHITE, LT_GRAY, GRAY, DK_GRAY, GRAY_BLACK, colorRGBA(0.0, 0.0, 0.5), colorRGBA(0.5, 0.0, 0.0)};
colorRGBA const trash_colors  [NUM_TRASH_COLORS  ] = {WHITE, WHITE, WHITE, WHITE, cream, vlt_yellow, LT_GRAY, colorRGBA(0.8, 0.6, 0.4)};
colorRGBA const pfloat_colors [NUM_PFLOAT_COLORS] = {WHITE, YELLOW, PINK, GREEN, ORANGE, LT_BLUE};
colorRGBA const LAMP_COLOR(1.0, 0.8, 0.6); // soft white
colorRGBA const WALL_LAMP_COLOR(1.0, 0.9, 0.8);
colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color; typical value to use
colorRGBA const DUCT_COLOR(WHITE);
colorRGBA const rat_color(GRAY); // make the rat's fur darker
colorRGBA const candle_color(0.95, 0.9, 0.75, 1.0); // cream
colorRGBA const GLASS_COLOR(0.8, 1.0, 0.9, 0.25);
colorRGBA const mall_tc_legs_color(BKGRAY);

unsigned const NUM_LOCK_COLORS = 8;
unsigned const MAX_LOCK_INDEX  = NUM_LOCK_COLORS + 2;
colorRGBA   const lock_colors     [NUM_LOCK_COLORS] = {WHITE, BLACK, RED, GREEN, BLUE, YELLOW, ORANGE, BROWN};
std::string const lock_color_names[NUM_LOCK_COLORS] = {"silver", "black", "red", "green", "blue", "yellow", "orange", "brown"};

inline colorRGBA gen_box_color(rand_gen_t &rgen) {return colorRGBA(rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0));} // add minor color variation
template<typename T> static bool check_vect_cube_contains_pt(vector<T> const &cubes, point const &pos);

class light_source;
class lmap_manager_t;
class building_nav_graph_t;
struct building_t;
class building_creator_t;
struct elevator_t;
struct escalator_t;
class brg_batch_draw_t;
typedef vector<vert_norm_comp_tc_color> vect_vnctcc_t;
struct sign_t;
struct city_flag_t;
struct door_t;
struct pipe_t;
struct tunnel_seg_t;
struct bldg_industrial_info_t;
typedef vector<point> vect_point;
typedef vector<sphere_t> vect_sphere_t;

struct bottle_params_t {
	std::string name, texture_fn;
	colorRGBA glass_color, liquid_color;
	float value, label_tscale;
	bottle_params_t(std::string const &n, std::string const &fn, colorRGBA const &gc, colorRGBA const &lc, float v, float ts) :
		name(n), texture_fn(fn), glass_color(gc), liquid_color(lc), value(v), label_tscale(ts) {}
};
enum {BOTTLE_TYPE_WATER=0, BOTTLE_TYPE_COKE, BOTTLE_TYPE_BEER, BOTTLE_TYPE_WINE, BOTTLE_TYPE_POISON, BOTTLE_TYPE_MEDS, NUM_BOTTLE_TYPES};
// Note: we could add colorRGBA(0.8, 0.9, 1.0, 0.4) for water bottles, but transparent objects require removing interior faces such as half of the sphere
bottle_params_t const bottle_params[NUM_BOTTLE_TYPES] = {
	bottle_params_t("bottle of water",    "interiors/arrowhead_logo.jpg", colorRGBA(0.4, 0.7, 1.0 ), colorRGBA(0.0,  0.5,  1.0, 0.0), 1.0,  1.0),
	bottle_params_t("bottle of Coke",     "interiors/coke_label.jpg",     colorRGBA(0.2, 0.1, 0.05), colorRGBA(0.22, 0.11, 0.06),     1.0,  1.0),
	bottle_params_t("bottle of beer",     "interiors/heineken_label.jpg", colorRGBA(0.1, 0.4, 0.1 ), colorRGBA(0.5,  0.4,  0.1 ),     3.0,  2.0),
	bottle_params_t("bottle of wine",     "interiors/wine_label.jpg",     BLACK,                     colorRGBA(0.4,  0.0,  0.1 ),     10.0, 2.0),
	bottle_params_t("bottle of poison",   "yuck.png",                     BLACK,                     BLACK,                           5.0,  2.0),
	bottle_params_t("bottle of medicine", "interiors/magenta_cross.png",  LT_BLUE,                   colorRGBA(0.2,  0.0,  0.7 ),     20.0, 1.0),
};

struct drink_can_params_t {
	std::string name, texture_fn;
	colorRGBA liquid_color;
	float value, tscale, tex_clip_y;
	drink_can_params_t(std::string const &n, std::string const &fn, colorRGBA const &lc, float v, float ts, float tcy) :
		name(n), texture_fn(fn), liquid_color(lc), value(v), tscale(ts), tex_clip_y(tcy) {}
};
enum {DRINK_CAN_TYPE_COKE=0, DRINK_CAN_TYPE_BEER, NUM_DRINK_CAN_TYPES};

drink_can_params_t const drink_can_params[NUM_DRINK_CAN_TYPES] = {
	drink_can_params_t("can of Coke", "interiors/coke_label.jpg",     colorRGBA(0.22, 0.11, 0.06), 1.0, 1.0, 0.16),
	drink_can_params_t("can of beer", "interiors/heineken_label.jpg", colorRGBA(0.5,  0.4,  0.1 ), 3.0, 2.0, 0.08),
};

struct ball_type_t {
	std::string name, tex_fname, nm_fname;
	float radius, density, value, weight, spec, shine, elastic, friction; // radius in inches, value in dollars, weight in pounds
	bool can_kick, hurts_zombie, breaks_glass;
	ball_type_t(std::string const &name_, std::string const &fn, std::string const &nm, float r, float d, float v, float w, bool ck, bool hz, bool bg,
		float spec_, float shine_, float e, float f) :
		name(name_), tex_fname(fn), nm_fname(nm), radius(r), density(d), value(v), weight(w), spec(spec_), shine(shine_),
		elastic(e), friction(f), can_kick(ck), hurts_zombie(hz), breaks_glass(bg) {}
};
enum {BALL_TYPE_SOCCER=0, BALL_TYPE_BASKET, BALL_TYPE_SOFT, BALL_TYPE_TENNIS, BALL_TYPE_BEACH, NUM_BALL_TYPES};

// name tex_fname nm_fname radius density value weight can_kick hurts_zombie breaks_glass spec shine elastic friction
ball_type_t const ball_types[NUM_BALL_TYPES] = {
	ball_type_t("soccer ball", "balls/soccer_ball_diffuse.png", "balls/soccer_ball_normal.png", 4.4, 0.50, 12.0, 0.90, 1, 1, 1, 0.4, 60.0, 1.0, 1.0),
	ball_type_t("basketball",  "balls/basketball.png",          "",                             4.7, 0.50, 15.0, 1.38, 1, 1, 1, 0.2, 40.0, 1.0, 1.0),
	ball_type_t("softball",    "balls/softball.jpg",            "",                             1.9, 1.20,  5.0, 0.40, 0, 1, 1, 0.1, 20.0, 0.8, 1.0), // balls/softball_bump.jpg    bad format
	ball_type_t("tennis ball", "balls/tennis_ball.jpg",         "",                             1.3, 0.75,  2.0, 0.13, 0, 1, 0, 0.0,  0.0, 1.0, 2.0), // balls/tennis_ball_bump.jpg bad format
	ball_type_t("beach ball",  "balls/beachball.jpg",           "",                            10.0, 0.01, 10.0, 0.10, 1, 0, 0, 0.5, 80.0, 0.8, 1.0)
};
ball_type_t const pool_ball_type("pool ball", "balls/pool_balls.png", "",                     1.125, 1.70,  2.0, 0.37, 0, 1, 1, 0.9, 100.0,0.5, 2.5);

class light_ix_assign_t {
	vector<pair<point2d<float>, unsigned>> cur;
	unsigned next_ix=0;
public:
	void next_room() {cur.clear();}
	unsigned get_next_ix() {return next_ix++;}
	unsigned get_ix_for_light(cube_t const &c, bool walls_not_shared=0);
};

struct building_occlusion_state_t {
	int exclude_bix=-1;
	bool skip_cont_camera=0;
	point pos;
	vector3d xlate;
	vector<cube_with_ix_t> building_ids;

	void init(point const &pos_, vector3d const &xlate_) {
		pos   = pos_;
		xlate = xlate_;
		building_ids.clear();
	}
};

class occlusion_checker_t {
	building_occlusion_state_t state;
public:
	vect_cube_t occluders;
	void set_exclude_bix(int exclude_bix) {state.exclude_bix = exclude_bix;}
	void set_exclude_camera_building() {state.skip_cont_camera = 1;}
	void set_camera(pos_dir_up const &pdu);
	bool is_occluded(cube_t const &c) const;
};

class occlusion_checker_noncity_t {
	building_occlusion_state_t state;
	building_creator_t const &bc;
public:
	bool query_is_for_light, extra_occluders_dim=0;
	vect_cube_t extra_occluders;

	occlusion_checker_noncity_t(building_creator_t const &bc_, bool for_light=0) : bc(bc_), query_is_for_light(for_light) {}
	void set_exclude_bix(int exclude_bix) {state.exclude_bix = exclude_bix;}
	void set_camera(pos_dir_up const &pdu, bool cur_building_only=0);
	bool is_occluded(cube_t const &c) const;
	vector3d const &get_xlate() const {return state.xlate;}
};

struct city_zone_t : public cube_t {
	float zval=0.0;
	bool is_park=0, is_residential=0;
	uint8_t street_dir=0; // encoded as 2*dim + dir + 1; 0 is unassigned
	unsigned nbuildings=0, capacity=0; // in number of buildings; 0 is unlimited
	unsigned max_floors=0; // 0=unlimited
	unsigned city_ix=0, street_num=0;
	int parent_plot_ix=-1; // if this is a sub-plot; -1 otherwise
	std::string address; // or just store road_name?

	city_zone_t() {}
	city_zone_t(cube_t const &c, float zval_, bool p, bool r, unsigned sdir, unsigned cap, int ppix, int cix, unsigned mf) :
		cube_t(c), zval(zval_), is_park(p), is_residential(r), street_dir(sdir), capacity(cap), max_floors(mf), city_ix(cix), parent_plot_ix(ppix) {}
	bool is_full() const {return (capacity > 0 && nbuildings >= capacity);}
};

typedef vector<city_zone_t> vect_city_zone_t;

struct tid_nm_pair_dstate_t {
	shader_t &s;
	int bmm_loc=-1;
	float bump_map_mag=1.0, crack_weight=0.0;
	bool no_set_texture=0;

	tid_nm_pair_dstate_t(shader_t &s_, bool no_set_texture_=0, float crack_weight_=0.0) : s(s_), crack_weight(crack_weight_), no_set_texture(no_set_texture_) {}
	void set_for_shader(float new_bump_map_mag);
	~tid_nm_pair_dstate_t();
};

struct tid_nm_pair_t { // size=32

	int tid=-1, nm_tid=-1; // Note: assumes each tid has only one nm_tid
	float tscale_x=1.0, tscale_y=1.0, txoff=0.0, tyoff=0.0, emissive=0.0;
	color_wrapper spec_color;
	unsigned char shininess=0; // Note: spec_mag is divided by 255.0
	bool shadowed   =0; // Note: doesn't directly affect rendering, only used for uniquing/operator==()
	bool shadow_only=0;
	bool transparent=0; // used to draw batched alpha blended materials last
	bool no_cracks  =0; // for basement crack effects

	tid_nm_pair_t() {}
	tid_nm_pair_t(int tid_, float txy=1.0, bool shadowed_=0, bool transparent_=0) : tid(tid_), nm_tid(FLAT_NMAP_TEX),
		tscale_x(txy), tscale_y(txy), shadowed(shadowed_), transparent(transparent_) {} // non-normal mapped 1:1 texture AR
	tid_nm_pair_t(int tid_, int nm_tid_, float tx, float ty, float xo=0.0, float yo=0.0, bool shadowed_=0, bool transparent_=0) :
		tid(tid_), nm_tid(nm_tid_), tscale_x(tx), tscale_y(ty), txoff(xo), tyoff(yo), shadowed(shadowed_), transparent(transparent_) {}
	void set_shininess(float shine) {shininess = (unsigned char)max(1, min(255, round_fp(shine)));}
	void set_specular(float mag, float shine) {set_specular_color(WHITE, mag, shine);}
	void set_specular_color(colorRGB const &color, float mag, float shine);
	void set_metal_specular(colorRGB const &color=WHITE) {set_specular_color(color, 0.8, 60.0);}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}

	bool is_compat_ignore_shadowed(tid_nm_pair_t const &t) const {
		return (tid == t.tid && nm_tid == t.nm_tid && emissive == t.emissive && shininess == t.shininess && transparent == t.transparent && spec_color == t.spec_color);
	}
	bool is_compatible(tid_nm_pair_t const &t) const {return (is_compat_ignore_shadowed(t) && shadowed == t.shadowed && shadow_only == t.shadow_only);}
	bool operator==(tid_nm_pair_t const &t) const {return (is_compatible(t) && tscale_x == t.tscale_x && tscale_y == t.tscale_y && txoff == t.txoff && tyoff == t.tyoff);}
	bool operator!=(tid_nm_pair_t const &t) const {return !operator==(t);}
	int get_nm_tid() const {return ((nm_tid < 0) ? FLAT_NMAP_TEX : nm_tid);}
	float get_drawn_tscale_x() const {return 2.0f*tscale_x;} // adjust for local vs. global space change
	float get_drawn_tscale_y() const {return 2.0f*tscale_y;} // adjust for local vs. global space change
	float get_emissive_val  () const;
	colorRGBA get_avg_color () const {return texture_color(tid);}
	tid_nm_pair_t get_scaled_version(float scale) const;
	static bool bind_reflection_shader();
	void   set_gl(tid_nm_pair_dstate_t &state) const;
	void unset_gl(tid_nm_pair_dstate_t &state) const;
	void toggle_transparent_windows_mode();
};

struct building_tex_params_t {
	tid_nm_pair_t side_tex, roof_tex; // exterior
	tid_nm_pair_t wall_tex, ceil_tex, floor_tex, house_ceil_tex, house_floor_tex, basement_floor_tex; // interior

	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0 || wall_tex.nm_tid >= 0 || ceil_tex.nm_tid >= 0 ||
		floor_tex.nm_tid >= 0 || house_ceil_tex.nm_tid >= 0 || house_floor_tex.nm_tid >= 0 || basement_floor_tex.nm_tid >= 0);}
};

struct color_range_t {
	float grayscale_rand=0.0;
	colorRGBA cmin=WHITE, cmax=WHITE; // alpha is unused?
	void gen_color(colorRGBA &color, rand_gen_t &rgen) const;
};

struct riser_pos_t : public sphere_t {
	bool has_hot=0, flow_dir=0, upper_floor=0, in_extb=0; // flow_dir: 0=out/down, 1=in/up
	riser_pos_t() {}
	riser_pos_t(point const &pos_, float radius_, bool hh=0, bool fd=0, bool uf=0, bool eb=0) : sphere_t(pos_, radius_), has_hot(hh), flow_dir(fd), upper_floor(uf), in_extb(eb) {}
};
typedef vector<riser_pos_t> vect_riser_pos_t;

struct building_mat_t : public building_tex_params_t {

	bool no_city=0, add_windows=0, add_wind_lights=0, no_walkways=0;
	unsigned min_levels=1, max_levels=1, min_sides=4, max_sides=4;
	float place_radius=0.0, max_delta_z=0.0, max_rot_angle=0.0, min_level_height=0.0, min_alt=-1000, max_alt=1000, house_prob=0.0, house_scale_min=1.0, house_scale_max=1.0;
	float split_prob=0.0, cube_prob=1.0, round_prob=0.0, asf_prob=0.0, min_fsa=0.0, max_fsa=0.0, min_asf=0.0, max_asf=0.0;
	float wind_xscale=1.0, wind_yscale=1.0, wind_xoff=0.0, wind_yoff=0.0;
	float floor_spacing=0.0, floorplan_wind_xscale=0.0; // these are derived values
	float apartment_prob=0.0;
	cube_t pos_range, prev_pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color; // exterior
	colorRGBA window_color=GRAY, wall_color=WHITE, ceil_color=WHITE, floor_color=LT_GRAY, house_ceil_color=WHITE, house_floor_color=WHITE;

	building_mat_t() : pos_range(-100,100,-100,100,0,0), sz_range(1,1,1,1,1,1) {}
	float gen_house_size_scale(rand_gen_t &rgen) const {return ((house_scale_min == house_scale_max) ? house_scale_min : rgen.rand_uniform(house_scale_min, house_scale_max));}
	void update_range(vector3d const &range_translate);
	void set_pos_range(cube_t const &new_pos_range) {prev_pos_range = pos_range; pos_range = new_pos_range;}
	void restore_prev_pos_range() {
		if (!prev_pos_range.is_all_zeros()) {pos_range = prev_pos_range;}
	}
	void finalize();
	float get_window_tx() const;
	float get_window_ty() const;
	float get_floor_spacing() const {return floor_spacing;}
	float get_floorplan_window_xscale() const {return floorplan_wind_xscale;}
};

struct building_params_t {

	bool flatten_mesh=0, has_normal_map=0, tex_mirror=0, tex_inv_y=0, tt_only=0, infinite_buildings=0, dome_roof=0, onion_roof=0;
	bool gen_building_interiors=1, add_city_interiors=0, enable_rotated_room_geom=0, add_secondary_buildings=0, add_office_basements=0, add_office_br_basements=0;
	bool put_doors_in_corners=0, cities_all_bldg_mats=0, small_city_buildings=0, add_door_handles=0, use_voronoise_cracks=0, add_basement_tunnels=0;
	unsigned num_place=0, num_tries=10, cur_prob=1, max_shadow_maps=32, buildings_rand_seed=0, max_ext_basement_hall_branches=4, max_ext_basement_room_depth=4;
	unsigned max_room_geom_gen_per_frame=1, max_office_basement_floors=2, max_mall_levels=2;
	float ao_factor=0.0, sec_extra_spacing=0.0, player_coll_radius_scale=1.0, interior_view_dist_scale=1.0;
	float window_width=0.0, window_height=0.0, window_xspace=0.0, window_yspace=0.0; // windows
	float wall_split_thresh=4.0, max_fp_wind_xscale=0.0, max_fp_wind_yscale=0.0, basement_water_level_min=0.0, basement_water_level_max=0.0; // interiors
	float open_door_prob=1.0, locked_door_prob=0.0, basement_prob_house=0.5, basement_prob_office=0.5, ball_prob=0.3, two_floor_retail_prob=0.0; // interior probabilities
	float split_stack_floorplan_prob=0.0, retail_floorplan_prob=0.0, mall_prob=0.5; // floorplan probabilities
	float glass_floor_alpha=GLASS_COLOR.A, max_altitude_all=0.0;
	// consistency probabilities of houses for cities and blocks
	float house_same_mat_prob =0.0, house_same_size_prob =0.0, house_same_geom_prob =0.0, house_same_per_city_prob =0.0;
	float office_same_mat_prob=0.0, office_same_size_prob=0.0, office_same_geom_prob=0.0, office_same_per_city_prob=0.0;
	// building people/AI params
	bool enable_people_ai=0, ai_target_player=1, ai_follow_player=0, allow_elevator_line=1, no_coll_enter_exit_elevator=1, show_player_model=0;
	unsigned ai_opens_doors=1; // 0=don't open doors, 1=only open if player closed door after path selection; 2=always open doors
	unsigned ai_player_vis_test=0; // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
	unsigned ai_sees_player_hide=2; // 0=doesn't see the player, 1=sees the player and waits outside the hiding spot, 2=opens the door and comes in
	unsigned people_per_office_min=0, people_per_office_max=0, people_per_house_min=0, people_per_house_max=0, elevator_capacity=1;
	unsigned player_model_ix=0;
	float ai_retreat_time=4.0, elevator_wait_time=60.0, use_elevator_prob=0.25, elevator_wait_recall_prob=0.5;
	float people_min_alpha=0.0;
	// building animal params
	unsigned num_rats_min=0, num_rats_max=0, min_attack_rats=0, num_spiders_min=0, num_spiders_max=0, num_snakes_min=0, num_snakes_max=0, num_insects_min=0, num_insects_max=0;
	float rat_speed   =0.0, rat_size_min   =0.5, rat_size_max   =1.0; // rats
	float spider_speed=0.0, spider_size_min=0.5, spider_size_max=1.0, spider_drawer_prob=0.0; // spiders
	float snake_speed =0.0, snake_size_min =0.5, snake_size_max =1.0; // snakes
	float insect_speed=0.0, insect_size_min=0.5, insect_size_max=1.0; // snakes
	// gameplay state
	float player_weight_limit=100.0;
	// materials
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix, mat_gen_ix_city, mat_gen_ix_nocity, mat_gen_ix_res; // {any, city_only, non_city, residential}
	vector<unsigned> rug_tids, picture_tids, desktop_tids, sheet_tids, paper_tids, food_box_tids, flag_tids, metal_tids;
	vector<std::string> food_box_names; // same size as food_box_tids
	// use for option reading
	int read_error=0;
	kw_to_val_map_t<bool     >  kwmb;
	kw_to_val_map_t<unsigned >  kwmu;
	kw_to_val_map_t<float    >  kwmf;
	kw_to_val_map_t<colorRGBA>  kwmc;
	kw_to_val_map_float_check_t kwmr;

	building_params_t(unsigned num=0) : num_place(num),
		kwmb(read_error, "buildings"), kwmu(read_error, "buildings"), kwmf(read_error, "buildings"), kwmc(read_error, "buildings"), kwmr(read_error, "buildings") {init_kw_maps();}
	bool parse_buildings_option(FILE *fp);
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
	bool building_people_enabled() const {return (people_per_office_max > 0 || people_per_house_max > 0);}
	bool windows_enabled  () const {return (window_width > 0.0 && window_height > 0.0 && window_xspace > 0.0 && window_yspace);} // all must be specified as nonzero
	bool gen_inf_buildings() const {return (infinite_buildings && world_mode == WMODE_INF_TERRAIN);}
	float get_window_width_fract () const {assert(windows_enabled()); return window_width /(window_width  + window_xspace);}
	float get_window_height_fract() const {assert(windows_enabled()); return window_height/(window_height + window_yspace);}
	float get_window_tx() const {assert(windows_enabled()); return 1.0f/(window_width  + window_xspace);}
	float get_window_ty() const {assert(windows_enabled()); return 1.0f/(window_height + window_yspace);}
	void add_cur_mat();
	void finalize();

	building_mat_t const &get_material(unsigned mat_ix) const {
		assert(mat_ix < materials.size());
		return materials[mat_ix];
	}
	vector<unsigned> const &get_mat_list(bool city_only, bool non_city_only, bool residential) const;
	unsigned choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only, bool residential) const;
	float get_max_house_size() const;
	void set_pos_range(cube_t const &pos_range);
	void restore_prev_pos_range();
private:
	void init_kw_maps();
	int read_building_texture(FILE *fp, std::string const &str, bool is_normal_map, int &error, bool check_filename=0, bool *no_cracks=nullptr);
	void read_texture_and_add_if_valid(FILE *fp, std::string const &str, int &error, vector<unsigned> &tids);
};

class building_draw_t;

struct building_geom_t { // describes the physical shape of a building
	unsigned num_sides;
	uint8_t door_sides[4]={}; // bit mask for 4 door sides, one per base part
	bool half_offset=0, is_pointed=0;
	float rot_sin=0.0, rot_cos=0.0, flat_side_amt=0.0, alt_step_factor=0.0, start_angle=0.0; // rotation in XY plane, around Z (up) axis
	//float roof_recess;

	building_geom_t(unsigned ns=4, float rs=0.0, float rc=1.0, bool pointed=0) : num_sides(ns), is_pointed(pointed), rot_sin(rs), rot_cos(rc) {}
	bool is_rotated() const {return (rot_sin != 0.0);}
	bool is_cube()    const {return (num_sides == 4);}
	bool use_cylinder_coll() const {return (num_sides > 8 && flat_side_amt == 0.0);} // use cylinder collision if not a cube, triangle, octagon, etc. (approximate)
	void do_xy_rotate    (point const &center, point &pos) const;
	void do_xy_rotate_inv(point const &center, point &pos) const;
	void do_xy_rotate_normal(point &n) const;
	void do_xy_rotate_normal_inv(point &n) const;
};

struct tquad_with_ix_t : public tquad_t {
	// roof {office, peak, hip, slope}, roof access cover, wall, house door, building front door, building back door, garage door, interior doors back face, office doors, roof doors
	enum {TYPE_ROOF_OFFICE=0, TYPE_ROOF_PEAK, TYPE_ROOF_HIP, TYPE_ROOF_SLOPE, TYPE_ROOF_ACC, TYPE_WALL, TYPE_HDOOR, TYPE_BDOOR, TYPE_BDOOR2, TYPE_GDOOR,
		TYPE_IDOOR, TYPE_IDOOR_IN, TYPE_ODOOR, TYPE_ODOOR_IN, TYPE_RDOOR, TYPE_RDOOR2, TYPE_RDOOR_IN, TYPE_HELIPAD, TYPE_SOLAR, TYPE_TRIM};
	bool is_roof         () const {return (type == TYPE_ROOF_OFFICE || type == TYPE_ROOF_PEAK || type == TYPE_ROOF_HIP || type == TYPE_ROOF_SLOPE);}
	bool is_building_door() const {return (type == TYPE_BDOOR || type == TYPE_BDOOR2);} // for office buildings
	bool is_exterior_door() const {return (type == TYPE_HDOOR || type == TYPE_GDOOR    || is_rooftop_door()  || is_building_door());}
	bool is_interior_door() const {return (type == TYPE_IDOOR || type == TYPE_IDOOR_IN || type == TYPE_ODOOR || type == TYPE_ODOOR_IN);}
	bool is_inside_face  () const {return (type == TYPE_IDOOR_IN || type == TYPE_ODOOR_IN || type == TYPE_RDOOR_IN);}
	bool is_rooftop_door () const {return (type == TYPE_RDOOR || type == TYPE_RDOOR2 || type == TYPE_RDOOR_IN);}

	unsigned type;
	tquad_with_ix_t(unsigned npts_=0, unsigned type_=TYPE_ROOF_PEAK) : tquad_t(npts_), type(type_) {}
	tquad_with_ix_t(tquad_t const &t, unsigned type_) : tquad_t(t), type(type_) {}
};
typedef vector<tquad_with_ix_t> vect_tquad_with_ix_t;

struct vertex_range_t { // size=12
	int draw_ix; // -1 is unset
	unsigned start, end;
	vertex_range_t() : draw_ix(-1), start(0), end(0) {}
	vertex_range_t(unsigned s, unsigned e, int ix=-1) : draw_ix(ix), start(s), end(e) {}
};

struct draw_range_t {
	// intended for building interiors, which don't have many materials; may need to increase MAX_DRAW_BLOCKS later; max used value appears to be 6 vrq and 2 vrt
	vertex_range_t vrq[MAX_DRAW_BLOCKS]; // quad     verts
	vertex_range_t vrt[MAX_DRAW_BLOCKS]; // triangle verts
};

// building types/functions; these are for primary buildings, not basements/rooms (such as malls or parking garages)
enum {BTYPE_UNSET=0, BTYPE_HOUSE, BTYPE_MULT_FAM, BTYPE_OFFICE, BTYPE_APARTMENT, BTYPE_HOTEL, BTYPE_HOSPITAL, BTYPE_PARKING, BTYPE_MALL, BTYPE_FACTORY,
	BTYPE_WAREHOUSE, BTYPE_POWERPLANT, BTYPE_SCHOOL, BTYPE_POLICE, BTYPE_FIRE_STAT, NUM_BUILDING_TYPES};
std::string const btype_names[NUM_BUILDING_TYPES] =
{"", "House", "Multi-Family House", "Office", "Apartments", "Hotel", "Hospital", "Parking", "Mall", "Factory", "Warehouse", "Power Plant", "School", "Police Station", "Fire Station"};
colorRGBA const  btype_colors[NUM_BUILDING_TYPES] =
{WHITE, WHITE, YELLOW,               WHITE,    GREEN,        GREEN,   BLUE,       BROWN,     ORANGE, RED,       RED,         RED,           PURPLE,   MAGENTA,          MAGENTA};
typedef uint8_t building_type_t;

enum { // room object types
	TYPE_NONE=0, TYPE_TABLE, TYPE_CHAIR, TYPE_STAIR, TYPE_STAIR_WALL, TYPE_ELEVATOR, TYPE_LIGHT, TYPE_RUG, TYPE_PICTURE, TYPE_WBOARD,
	TYPE_BOOK, TYPE_BCASE, TYPE_TCAN, TYPE_DESK, TYPE_BED, TYPE_WINDOW, TYPE_BLOCKER, TYPE_COLLIDER, TYPE_CUBICLE, TYPE_STALL,
	TYPE_SIGN, TYPE_COUNTER, TYPE_CABINET, TYPE_KSINK, TYPE_BRSINK, TYPE_PLANT, TYPE_DRESSER, TYPE_NIGHTSTAND, TYPE_FLOORING, TYPE_CLOSET,
	TYPE_WALL_TRIM, TYPE_RAILING, TYPE_CRATE, TYPE_BOX, TYPE_MIRROR, TYPE_SHELVES, TYPE_KEYBOARD, TYPE_SHOWER, TYPE_RDESK, TYPE_BOTTLE,
	TYPE_WINE_RACK, TYPE_COMPUTER, TYPE_MWAVE, TYPE_PAPER, TYPE_BLINDS, TYPE_PEN, TYPE_PENCIL, TYPE_PAINTCAN, TYPE_LG_BALL, TYPE_HANGER_ROD,
	TYPE_DRAIN, TYPE_MONEY, TYPE_PHONE, TYPE_TPROLL, TYPE_SPRAYCAN, TYPE_MARKER, TYPE_BUTTON, TYPE_CRACK, TYPE_SWITCH, TYPE_PLATE,
	TYPE_LAPTOP, TYPE_FPLACE, TYPE_LBASKET, TYPE_WHEATER, TYPE_TAPE, TYPE_OUTLET, TYPE_PG_WALL, TYPE_PG_PILLAR, TYPE_PG_BEAM, TYPE_PARK_SPACE,
	TYPE_RAMP, TYPE_PIPE, TYPE_CURB, TYPE_BRK_PANEL, TYPE_VENT, TYPE_BREAKER, TYPE_FURNACE, TYPE_ATTIC_DOOR, TYPE_CHIMNEY, TYPE_DUCT,
	TYPE_TOY, TYPE_DRESS_MIR, TYPE_PAN, TYPE_VASE, TYPE_URN, TYPE_FCABINET, TYPE_STAPLER, TYPE_WIND_SILL, TYPE_BALCONY, TYPE_SPRINKLER,
	TYPE_FEXT_MOUNT, TYPE_FEXT_SIGN, TYPE_PIZZA_BOX, TYPE_PIZZA_TOP, TYPE_TEESHIRT, TYPE_PANTS, TYPE_BLANKET, TYPE_SERVER, TYPE_EXT_STEP, TYPE_DBG_SHAPE,
	TYPE_POOL_BALL, TYPE_POOL_CUE, TYPE_WALL_MOUNT, TYPE_POOL_TILE, TYPE_POOL_FLOAT, TYPE_BENCH, TYPE_DIV_BOARD, TYPE_FALSE_DOOR, TYPE_FLASHLIGHT, TYPE_CANDLE,
	TYPE_CAMERA, TYPE_CLOCK, TYPE_DOWNSPOUT, TYPE_SHELFRACK, TYPE_CHIM_CAP, TYPE_FOOD_BOX, TYPE_SAFE, TYPE_LADDER, TYPE_CHECKOUT, TYPE_FISHTANK,
	TYPE_LAVALAMP, TYPE_SHOWERTUB, TYPE_TRASH, TYPE_VALVE, TYPE_METAL_BAR, TYPE_OFF_PILLAR, TYPE_DRINK_CAN, TYPE_CONF_TABLE, TYPE_INT_WINDOW, TYPE_INT_LADDER,
	TYPE_MACHINE, TYPE_BUCKET, TYPE_SPIWEB, TYPE_TREE, TYPE_THEFT_SENS, TYPE_ELEC_WIRE, TYPE_ERASER, TYPE_DWASHER, TYPE_PET_CAGE, TYPE_IBEAM,
	TYPE_CATWALK, TYPE_VANITY, TYPE_CHEM_TANK, TYPE_HVAC_UNIT, TYPE_WARN_LIGHT, TYPE_GAUGE, TYPE_PALLET, TYPE_SHELF_WALL, TYPE_VENDING, TYPE_MED_CAB,
	TYPE_LOCKER,
	/* these next ones are all 3D models - see logic in room_object_t::is_obj_model_type() */
	TYPE_TOILET, TYPE_SINK, TYPE_TUB, TYPE_FRIDGE, TYPE_STOVE, TYPE_TV, TYPE_MONITOR, TYPE_COUCH, TYPE_OFF_CHAIR, TYPE_URINAL,
	TYPE_LAMP, TYPE_WASHER, TYPE_DRYER, TYPE_KEY, TYPE_HANGER, TYPE_CLOTHES, TYPE_FESCAPE, TYPE_WALL_LAMP, TYPE_CUP, TYPE_TOASTER,
	TYPE_HOOD, TYPE_RCHAIR, TYPE_SILVER, TYPE_TOY_MODEL, TYPE_CEIL_FAN, TYPE_FIRE_EXT, TYPE_FOLD_SHIRT, TYPE_PLANT_MODEL, TYPE_POOL_TABLE, TYPE_POOL_LAD,
	TYPE_BAR_STOOL, TYPE_PADLOCK, TYPE_CASHREG, TYPE_WFOUNTAIN, TYPE_BANANA, TYPE_BAN_PEEL, TYPE_CONF_PHONE, TYPE_SHOE, TYPE_SHOEBOX, TYPE_VENT_FAN,
	TYPE_HOSP_BED, TYPE_HOSP_CURT, TYPE_FORKLIFT, TYPE_WHEELCHAIR, TYPE_OP_TABLE, TYPE_TROLLEY,
	/* shared with city objects */
	TYPE_GBIKE, TYPE_XFORMER, TYPE_US_FLAG, TYPE_BLDG_FOUNT,
	/* animals; bird is only used for pet stores */
	TYPE_RAT, TYPE_ROACH, TYPE_SPIDER, TYPE_SNAKE, TYPE_INSECT, TYPE_FISH, TYPE_BIRD,
	NUM_ROBJ_TYPES};
typedef uint8_t room_object;

// room object and stairs shapes
enum {SHAPE_CUBE=0, SHAPE_CYLIN, SHAPE_SPHERE, SHAPE_STAIRS_U, SHAPE_STAIRS_L, SHAPE_STAIRS_FAN, SHAPE_TALL, SHAPE_SHORT, SHAPE_ANGLED, SHAPE_VERT_TORUS};
typedef uint8_t room_obj_shape;

// room types
enum {RTYPE_NOTSET=0, RTYPE_HALL, RTYPE_STAIRS, RTYPE_OFFICE, RTYPE_BATH, RTYPE_MENS, RTYPE_WOMENS, RTYPE_BED, RTYPE_KITCHEN, RTYPE_LIVING,
	  RTYPE_DINING, RTYPE_STUDY, RTYPE_ENTRY, RTYPE_LIBRARY, RTYPE_STORAGE, RTYPE_GARAGE, RTYPE_SHED, RTYPE_LOBBY, RTYPE_LAUNDRY, RTYPE_CARD,
	  RTYPE_PLAY, RTYPE_ART, RTYPE_UTILITY, RTYPE_PARKING, RTYPE_RAMP_EXIT, RTYPE_ATTIC, RTYPE_MASTER_BED, RTYPE_UNFINISHED, RTYPE_SERVER, RTYPE_POOL,
	  RTYPE_SWIM, RTYPE_SECURITY, RTYPE_LOUNGE, RTYPE_COMMON, RTYPE_BACKROOMS, RTYPE_RETAIL, RTYPE_ELEVATOR, RTYPE_CONF, RTYPE_MACHINE, RTYPE_INTERR,
	  RTYPE_ELEV_EQUIP, RTYPE_STORE, RTYPE_MALL, RTYPE_RESTAURANT, RTYPE_FACTORY, RTYPE_WAREHOUSE, RTYPE_HOS_BED, RTYPE_HOS_OR, RTYPE_HOS_EXAM, RTYPE_CLASS,
	  RTYPE_WAITING, NUM_RTYPES};
typedef uint8_t room_type;

inline bool is_bathroom (room_type   const rtype) {return (rtype == RTYPE_BATH || rtype == RTYPE_MENS || rtype == RTYPE_WOMENS);}
inline bool is_ball_type(room_object const type ) {return (type == TYPE_LG_BALL || type == TYPE_POOL_BALL);}

// full room names for UI display
std::string const room_names[NUM_RTYPES] =
	{"Unset", "Hallway", "Stairs", "Office", "Bathroom", "Men's Restroom", "Women's Restroom", "Bedroom", "Kitchen", "Living Room",
	 "Dining Room", "Study", "Entryway", "Library", "Storage Room", "Garage", "Shed", "Lobby", "Laundry Room", "Card Room",
	 "Play Room", "Art Room", "Utility Room", "Parking Garage", "Ramp Exit", "Attic", "Master Bedroom", "Unfinished Room", "Server Room", "Pool Room",
	 "Swimming Pool Room", "Security Room", "Lounge", "Common Room", "Backrooms", "Retail", "Elevator", "Conference Room", "Machine Room", "Interrogation Room",
	 "Elev Equip Room", "Store", "Mall Concourse", "Restaurant", "Factory Floor", "Warehouse", "Hospital Bedroom", "Operating Room", "Exam Room", "Classroom",
	 "Waiting Room"
};
// short room names for elevator buttons (should be <= 8 characters)
std::string const room_names_short[NUM_RTYPES] =
	{"", "Hall", "Stairs", "Office", "Bath", "Men", "Women", "Bed", "Kitchen", "Living",
	"Dining", "Study", "Entry", "Library", "Storage", "Garage", "Shed", "Lobby", "Laundry", "Card",
	"Play", "Art", "Utility", "Garage", "Ramp", "Attic", "Bed", "", "Server", "Pool",
	"Swim", "Security", "Lounge", "Common", "Basement", "Retail", "Elevator", "Conference", "Machine", "Dungeon",
	"Equipment", "Store", "Mall", "Restaurant", "Factory", "Warehouse", "Bedroom", "OR", "Exam", "Class",
	"Waiting"
};

unsigned const room_priorities[NUM_RTYPES] = { // for breaker labels
	0, 2, 1, 1, 2, 2, 2, 3, 3, 3,
	3, 1, 3, 3, 2, 3, 3, 2, 2, 2,
	2, 2, 3, 3, 0, 3, 3, 0, 4, 3,
	4, 4, 4, 0, 1, 2, 1, 3, 2, 2,
	1, 3, 4, 3, 1, 1, 2, 3, 3, 2,
	2
};

// store types, for use with object placement and naming
enum {STORE_OTHER=0, STORE_CLOTHING, STORE_FOOD, STORE_BOOK, STORE_RETAIL, STORE_FURNITURE, STORE_PETS, STORE_APPLIANCE, STORE_SHOE, NUM_STORE_TYPES};
enum {RETAIL_BOXED=0, RETAIL_FOOD, RETAIL_HOUSE_GOODS, RETAIL_KITCHEN, RETAIL_ELECTRONICS, NUM_RETAIL_CAT};
std::string const store_type_strs [NUM_STORE_TYPES] = {"", "clothing", "food", "book", "retail", "furniture", "pet", "appliance", "shoe"};
std::string const srack_categories[NUM_RETAIL_CAT ] = {"boxed items", "food", "household goods", "kitchen", "electronics"};

enum {SHAPE_STRAIGHT=0, SHAPE_U, SHAPE_WALLED, SHAPE_WALLED_SIDES, SHAPE_RAMP, SHAPE_L, SHAPE_FAN}; // stairs shapes; SHAPE_FAN is used for mall entrances
typedef uint8_t stairs_shape;

enum {ROOM_WALL_INT=0, ROOM_WALL_SEP, ROOM_WALL_EXT, ROOM_WALL_BASEMENT};
enum {FLOORING_MARBLE=0, FLOORING_TILE, FLOORING_CONCRETE, FLOORING_CARPET, FLOORING_WOOD}; // Note: not all are used
enum {MAT_TYPE_STATIC=0, MAT_TYPE_SMALL, MAT_TYPE_DYNAMIC, MAT_TYPE_DETAIL, MAT_TYPE_DOORS, MAT_TYPE_LIGHTS, MAT_TYPE_TEXT}; // building_room_geom_t material types; max is 8
enum {FTYPE_NONE=0, FTYPE_BASEMENT, FTYPE_ATTIC}; // for furnace
enum {ATTIC_TYPE_RAFTERS=0, ATTIC_TYPE_FIBERGLASS, ATTIC_TYPE_WOOD, ATTIC_TYPE_PLASTER, NUM_ATTIC_TYPES};
enum {BLDG_COLL_NONE=0, BLDG_COLL_SIDE, BLDG_COLL_ROOF, BLDG_COLL_DETAIL, BLDG_COLL_DRIVEWAY, BLDG_COLL_FENCE, BLDG_COLL_SKYLIGHT};
enum {PLACED_TOILET=1, PLACED_SINK=2, PLACED_TUB=4, PLACED_SHOWER=8}; // for bathroom objects

// Note: when adding an entry here, must also add a string to model_opt_names in city_building_params.cpp
enum {/*building models*/ OBJ_MODEL_TOILET=0, OBJ_MODEL_SINK, OBJ_MODEL_TUB, OBJ_MODEL_FRIDGE, OBJ_MODEL_STOVE, OBJ_MODEL_TV, OBJ_MODEL_MONITOR/*unused*/, OBJ_MODEL_COUCH,
	OBJ_MODEL_OFFICE_CHAIR, OBJ_MODEL_URINAL, OBJ_MODEL_LAMP, OBJ_MODEL_WASHER, OBJ_MODEL_DRYER, OBJ_MODEL_KEY, OBJ_MODEL_HANGER, OBJ_MODEL_CLOTHES, OBJ_MODEL_FESCAPE,
	OBJ_MODEL_WALL_LAMP, OBJ_MODEL_CUP, OBJ_MODEL_TOASTER, OBJ_MODEL_HOOD, OBJ_MODEL_RCHAIR, OBJ_MODEL_SILVER, OBJ_MODEL_TOY, OBJ_MODEL_CEIL_FAN, OBJ_MODEL_FIRE_EXT,
	OBJ_MODEL_FOLD_SHIRT, OBJ_MODEL_PLANT, OBJ_MODEL_POOL_TABLE, OBJ_MODEL_POOL_LAD, OBJ_MODEL_BAR_STOOL, OBJ_MODEL_PADLOCK, OBJ_MODEL_CASHREG, OBJ_MODEL_WFOUNTAIN,
	OBJ_MODEL_BANANA, OBJ_MODEL_BAN_PEEL, OBJ_MODEL_PHONE, OBJ_MODEL_SHOE, OBJ_MODEL_SHOEBOX, OBJ_MODEL_VENT_FAN, OBJ_MODEL_HOSP_BED, OBJ_MODEL_HOSP_CURT, OBJ_MODEL_FORKLIFT,
	OBJ_MODEL_WHEELCHAIR, OBJ_MODEL_OP_TABLE, OBJ_MODEL_TROLLEY,
	OBJ_MODEL_GBIKE/*unused*/, OBJ_MODEL_XFMR/*unused*/, OBJ_MODEL_US_FLAG/*unused*/, OBJ_MODEL_BLDG_FOUNT/*unused*/,
	/*animal models*/ OBJ_MODEL_RAT, OBJ_MODEL_ROACH,
	/*building non-room objects*/ OBJ_MODEL_DOOR_HANDLE,
	/*city models*/ OBJ_MODEL_FHYDRANT, OBJ_MODEL_SUBSTATION, OBJ_MODEL_MAILBOX, OBJ_MODEL_UMBRELLA, OBJ_MODEL_PIGEON, OBJ_MODEL_FOUNTAIN, OBJ_MODEL_BIRD_ANIM, OBJ_MODEL_FLAG,
	OBJ_MODEL_BICYCLE, OBJ_MODEL_SWINGSET, OBJ_MODEL_TRAMPOLINE, OBJ_MODEL_DUMPSTER, OBJ_MODEL_BIG_UMBRELLA, OBJ_MODEL_FLOWER, OBJ_MODEL_DECK_CHAIR, OBJ_MODEL_PICNIC,
	OBJ_MODEL_WIND_TUR, NUM_OBJ_MODELS};

enum {PART_EFFECT_NONE=0, PART_EFFECT_SPARK, PART_EFFECT_CLOUD, PART_EFFECT_SMOKE, PART_EFFECT_SPLASH, PART_EFFECT_BUBBLE, PART_EFFECT_DROPLET, NUM_PART_EFFECTS};
enum {PIPE_TYPE_SEWER=0, PIPE_TYPE_CW, PIPE_TYPE_HW, PIPE_TYPE_GAS, NUM_PIPE_TYPES};

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};

// object flags
unsigned const RO_FLAG_LIT     = 0x01; // light is on
unsigned const RO_FLAG_TOS     = 0x02; // at top of stairs; used for railings and lights
unsigned const RO_FLAG_IN_WH   = 0x02; // in warehouse; aliased with RO_FLAG_TOS
unsigned const RO_FLAG_RSTAIRS = 0x04; // in a room with stairs
unsigned const RO_FLAG_INVIS   = 0x08; // invisible
unsigned const RO_FLAG_NOCOLL  = 0x10; // no collision detection
unsigned const RO_FLAG_OPEN    = 0x20; // open, for elevators, closet doors, bathroom stalls, phones, etc.
unsigned const RO_FLAG_NODYNAM = 0x40; // for light shadow maps
unsigned const RO_FLAG_INTERIOR= 0x80; // applies to containing room
// object flags, second byte
unsigned const RO_FLAG_EMISSIVE= 0x0100; // for signs, lights, and phones
unsigned const RO_FLAG_HANGING = 0x0200; // for signs, blinds, hangers, shirts, beams, walls, and balconies; treated as "folding" for closet doors
unsigned const RO_FLAG_ADJ_LO  = 0x0400; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_HI  = 0x0800; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_BOT = 0x1000; // for door trim/railings/ext steps/etc.
unsigned const RO_FLAG_ADJ_TOP = 0x2000; // for door trim/railings
unsigned const RO_FLAG_IS_HOUSE= 0x4000; // used for mirror reflections, shelves, tables, desks, beds, closets, and false doors
unsigned const RO_FLAG_RAND_ROT= 0x8000; // random rotation; used for office chairs, papers, pictures, cups, and balls
unsigned const RO_FLAG_UNTEXTURED= 0x1000; // for shirts, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_FROM_SET  = 0x1000; // for books,  aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_HAS_VOL_IX= 0x2000; // for books,  aliased with RO_FLAG_ADJ_TOP
unsigned const RO_FLAG_FOR_CAR   = 0x1000; // for car blockers, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_WALKWAY   = 0x1000; // for walkway objects (outside of buildings), aliased with RO_FLAG_ADJ_BOT
// object flags, third byte, for pickup/interact state
unsigned const RO_FLAG_IN_HALLWAY= 0x010000; // for attic doors and trashcans
unsigned const RO_FLAG_IN_MALL   = 0x010000; // for mall chairs, tables, benches, trashcans, etc.; aliased with RO_FLAG_IN_HALLWAY and RO_FLAG_IN_FACTORY
unsigned const RO_FLAG_IN_FACTORY= 0x010000; // for machines, ladders, etc.; aliased with RO_FLAG_IN_HALLWAY and RO_FLAG_IN_MALL
unsigned const RO_FLAG_IN_ATTIC  = 0x020000; // in attic
unsigned const RO_FLAG_HAS_EXTRA = 0x040000; // used for counter backsplash, ext wall trim, desks with comp monitors and kbds, books on glass tables, hotel closets, paper towels
unsigned const RO_FLAG_EXTERIOR  = 0x080000; // for signs, window trim, etc.
unsigned const RO_FLAG_EXPANDED  = 0x100000; // for shelves, closets, boxes, and mirrors
unsigned const RO_FLAG_WAS_EXP   = 0x200000; // for objects in/on shelves, closets, drawers, cabinets, shelfracks, and books
unsigned const RO_FLAG_ROTATING  = 0x400000; // for office chairs and clothes on hangers
unsigned const RO_FLAG_IN_CLOSET = 0x800000; // for closet lights and light switches
unsigned const RO_FLAG_ON_SRACK  = 0x800000; // on shelf rack; aliased with RO_FLAG_IN_CLOSET/RO_FLAG_ON_FLOOR
unsigned const RO_FLAG_NONEMPTY  = 0x040000; // for microwaves, shelves, shelfracks, and cups, aliased with RO_FLAG_HAS_EXTRA
unsigned const RO_FLAG_ON_FLOOR  = 0x800000; // for books, fallen objects, upper floor shelfracks, etc., aliased with RO_FLAG_IN_CLOSET/RO_FLAG_ON_SRACK
unsigned const RO_FLAG_BROKEN2   = 0x040000; // for lights that are completely broken, aliased with RO_FLAG_HAS_EXTRA and RO_FLAG_NONEMPTY
unsigned const RO_FLAG_PLCOLL    = 0x040000; // player collidable, for chairs, aliased with RO_FLAG_HAS_EXTRA
// object flags, fourth byte
unsigned const RO_FLAG_DYNAMIC  = 0x01000000; // dynamic object (balls, elevators, etc.)
unsigned const RO_FLAG_DSTATE   = 0x02000000; // this object has dynamic state
unsigned const RO_FLAG_NO_CONS  = 0x04000000; // this object is not consumable (bottles)
unsigned const RO_FLAG_NO_POWER = 0x04000000; // unpowered; related to circuit breakers, aliased with RO_FLAG_NO_CONS
unsigned const RO_FLAG_IS_ACTIVE= 0x08000000; // active, for sinks, tubs, buttons, pool balls, shower curtains, etc.
unsigned const RO_FLAG_USED     = 0x10000000; // used by the player (spraypaint, marker, etc.); used by parking spaces to indicate cars
unsigned const RO_FLAG_IN_ELEV  = 0x20000000; // for elevator lights, buttons, and flooring; aliased with RO_FLAG_BACKROOM and RO_FLAG_IN_POOL
unsigned const RO_FLAG_BACKROOM = 0x20000000; // in backrooms, for walls and pillars; aliased with RO_FLAG_IN_ELEV and RO_FLAG_IN_POOL
unsigned const RO_FLAG_IN_POOL  = 0x20000000; // for stairs, railings, and drains; aliased with RO_FLAG_IN_ELEV and RO_FLAG_BACKROOM
unsigned const RO_FLAG_BROKEN   = 0x40000000; // for TVs, monitors, flickering lights, and ond computers; maybe can use for windows
unsigned const RO_FLAG_MOVED    = 0x80000000; // for player push/pull

// reflection pass flags
unsigned const REF_PASS_ENABLED     = 0x01;
unsigned const REF_PASS_INTERIOR    = 0x02;
unsigned const REF_PASS_HOUSE       = 0x04;
unsigned const REF_PASS_WATER       = 0x08;
unsigned const REF_PASS_EXTB        = 0x10;
unsigned const REF_PASS_NO_MIRROR   = 0x20;
unsigned const REF_PASS_INT_ONLY    = 0x40;
unsigned const REF_PASS_GLASS_FLOOR = 0x80;

inline vector3d vector_from_dim_dir(int dim, bool dir) {
	vector3d v;
	v[dim] = (dir ? 1.0 : -1.0);
	return v;
}

struct bldg_obj_type_t {
	bool player_coll=0, ai_coll=0, rat_coll=0, pickup=0, attached=0, is_model=0;
	uint8_t lg_sm=0; // 0=neither (model or detail), 1=large item, 2=small item, 3=split into large and small
	float value=0.0, weight=0.0;
	unsigned capacity=0; // for consumable/usable objects
	std::string name;

	bldg_obj_type_t() {}
	bldg_obj_type_t(bool pc, bool ac, bool rc, bool pu, bool at, bool im, uint8_t ls, float v, float w, std::string const &n, unsigned cap=0) :
		player_coll(pc), ai_coll(ac), rat_coll(rc), pickup(pu), attached(at), is_model(im), lg_sm(ls), value(v), weight(w), capacity(cap), name(n) {}
	void update_modified_flags_for_type(bldg_obj_type_t t) {lg_sm |= t.lg_sm; is_model |= t.is_model; ai_coll |= t.ai_coll;}
};

struct oriented_cube_t : public cube_t {
	bool dim=0, dir=0;
	oriented_cube_t() {}
	oriented_cube_t(cube_t const &c, bool dim_, bool dir_) : cube_t(c), dim(dim_), dir(dir_) {}
	float get_length() const {return get_sz_dim( dim);}
	float get_width () const {return get_sz_dim(!dim);}
	float get_height() const {return dz();}
};

uint16_t const sink_water_state_bit(1 << 15); // MSB of uint16_t; hopefully doesn't conflict with doors or drawers

struct room_object_t : public oriented_cube_t { // size=68
	uint8_t taken_level=0;
	// Note: state_flags is used for drawer was_opened state, railing num_stairs, pool balls, and sink/tub/shower water
	uint16_t room_id=0, obj_id=0, drawer_flags=0, item_flags=0, state_flags=0;
	room_object type=TYPE_NONE; // 8-bit
	room_obj_shape shape=SHAPE_CUBE; // 8-bit
	unsigned flags=0;
	float light_amt=1.0;
	colorRGBA color;

	room_object_t() {}
	room_object_t(cube_t const &c, room_object type_, uint16_t rid, bool dim_=0, bool dir_=0, unsigned f=0, float light=1.0,
		room_obj_shape shape_=SHAPE_CUBE, colorRGBA const color_=WHITE, uint16_t iflags=0) :
		oriented_cube_t(c, dim_, dir_), room_id(rid), item_flags(iflags), type(type_), shape(shape_), flags(f), light_amt(light), color(color_)
	{check_normalized();}
	void check_normalized() const;
	// treat {drawer_flags, item_flags, state_flags} as a single 48-bit flags; the upper 16 bits of set_combined_flags() argument are ignored
	uint64_t get_combined_flags() const {return (uint64_t)drawer_flags + ((uint64_t)item_flags << 16) + ((uint64_t)state_flags << 32);}
	void set_combined_flags(uint64_t v) {drawer_flags = (v & 0xFFFF); item_flags = ((v >> 16) & 0xFFFF); state_flags = ((v >> 32) & 0xFFFF);}
	static uint64_t get_book_ix_mask(unsigned ix) {return (uint64_t)1 << (ix % MAX_BCASE_BOOKS);} // up to 48 books
	bool is_valid   () const {return  (type != TYPE_NONE);}
	bool is_lit     () const {return  (flags & RO_FLAG_LIT);}
	bool is_powered () const {return !(flags & RO_FLAG_NO_POWER);}
	bool is_light_on() const {return  (is_lit() && is_powered());}
	bool has_stairs () const {return  (flags & RO_FLAG_RSTAIRS);}
	bool is_visible () const {return !(flags & RO_FLAG_INVIS);}
	bool no_coll    () const {return  (flags & RO_FLAG_NOCOLL);}
	bool is_interior() const {return  (flags & RO_FLAG_INTERIOR);}
	bool is_open    () const {return  (flags & RO_FLAG_OPEN);}
	bool is_house   () const {return  (flags & RO_FLAG_IS_HOUSE);}
	bool was_expanded() const{return  (flags & RO_FLAG_WAS_EXP);}
	bool obj_expanded() const{return  (flags & RO_FLAG_EXPANDED);}
	bool is_dynamic () const {return  (flags & RO_FLAG_DYNAMIC);}
	bool has_dstate () const {return  (flags & RO_FLAG_DSTATE);}
	bool is_moving  () const {return (is_dynamic() && has_dstate());}
	bool was_moved  () const {return  (flags & RO_FLAG_MOVED);}
	bool is_broken  () const {return  (flags & RO_FLAG_BROKEN );}
	bool is_broken2 () const {return  (flags & RO_FLAG_BROKEN2);} // fully broken
	bool is_active  () const {return  (flags & RO_FLAG_IS_ACTIVE);}
	bool is_used    () const {return  (flags & RO_FLAG_USED);}
	bool is_hanging () const {return  (flags & RO_FLAG_HANGING);}
	bool in_elevator() const {return  (flags & RO_FLAG_IN_ELEV);}
	bool in_closet  () const {return  (flags & RO_FLAG_IN_CLOSET);}
	bool in_attic   () const {return  (flags & RO_FLAG_IN_ATTIC);}
	bool in_mall    () const {return  (flags & RO_FLAG_IN_MALL);}
	bool in_factory () const {return  (flags & RO_FLAG_IN_FACTORY);}
	bool in_warehouse()const {return  (flags & RO_FLAG_IN_WH);}
	bool in_hallway () const {return  (flags & RO_FLAG_IN_HALLWAY);}
	bool is_exterior() const {return  (flags & RO_FLAG_EXTERIOR);}
	bool rotates    () const {return  (flags & RO_FLAG_RAND_ROT);}
	bool is_on_floor() const {return  (flags & RO_FLAG_ON_FLOOR);}
	bool is_nonempty() const {return  (flags & RO_FLAG_NONEMPTY);}
	bool has_lid    () const {return  (flags & RO_FLAG_ADJ_TOP);} // for fishtanks
	bool has_extra  () const {return  (flags & RO_FLAG_HAS_EXTRA);}
	bool on_warehouse_floor() const {return (in_warehouse() && (flags & RO_FLAG_ON_FLOOR));}
	bool is_crate_or_box() const {return (type == TYPE_CRATE || type == TYPE_BOX);}
	bool is_on_srack() const {return (was_expanded() && (flags & RO_FLAG_ON_SRACK));} // Note: RO_FLAG_ON_SRACK is aliased with other flags, so also check was_expanded
	bool is_floor_clutter() const {return ((is_a_drink() || type == TYPE_TRASH) && is_on_floor());}
	bool is_light_type() const {return (type == TYPE_LIGHT || (type == TYPE_LAMP && !was_expanded() && !in_attic()));} // light, or lamp not in closet
	bool is_sink_type () const {return (type == TYPE_SINK || type == TYPE_KSINK || type == TYPE_BRSINK || type == TYPE_VANITY);}
	bool is_obj_model_type() const {return (type >= TYPE_TOILET && type < NUM_ROBJ_TYPES);}
	bool is_small_closet() const {return (type == TYPE_CLOSET && get_width() < 1.2*dz());}
	bool is_bottle_empty() const {return ((obj_id & 192) == 192);} // empty if both bits 6 and 7 are set; also applies to drink cans
	bool desk_has_drawers()const {return (((room_id & 3) || has_extra()) && get_width() > 1.5*get_depth());} // 75% of the time, if wide enough
	bool is_glass_table () const {return (type == TYPE_TABLE && is_house() && (obj_id & 1));} // 50% chance if in a house
	bool is_parked_car  () const {return (type == TYPE_COLLIDER && (flags & RO_FLAG_FOR_CAR));}
	bool is_sloped_ramp () const {return (type == TYPE_RAMP || (type == TYPE_POOL_TILE && shape == SHAPE_ANGLED));}
	bool light_is_out   () const {return ((is_broken() || is_broken2()) && !is_open());} // only makes sense to call on lights
	bool is_mirror      () const {return (type == TYPE_MIRROR || type == TYPE_DRESS_MIR || (type == TYPE_MED_CAB && !has_extra()));}
	bool is_tv_or_monitor() const {return (type == TYPE_TV || type == TYPE_MONITOR);}
	bool is_player_collidable() const;
	bool can_use        () const;
	bool is_interactive () const {return (has_dstate() || can_use());}
	bool is_medicine    () const {return (type == TYPE_BOTTLE && get_bottle_type() == BOTTLE_TYPE_MEDS && !is_bottle_empty());}
	bool can_place_onto () const;
	bool is_floor_collidable () const;
	bool is_spider_collidable() const;
	bool is_collidable(bool for_spider) const {return (for_spider ? is_spider_collidable() : is_floor_collidable());}
	bool is_vert_cylinder() const;
	bool is_round  () const {return (shape == SHAPE_CYLIN || shape == SHAPE_SPHERE || shape == SHAPE_VERT_TORUS);}
	bool is_a_drink() const {return (type == TYPE_BOTTLE || type == TYPE_DRINK_CAN);}
	bool is_pet_container() const {return (type == TYPE_FISHTANK || type == TYPE_PET_CAGE);}
	unsigned get_bottle_type   () const {return ((obj_id&63) % NUM_BOTTLE_TYPES   );} // first 6 bits are bottle type
	unsigned get_drink_can_type() const {return ((obj_id&63) % NUM_DRINK_CAN_TYPES);} // first 6 bits are drink can type
	unsigned get_orient() const {return (2*dim + dir);}
	unsigned get_num_shelves() const;
	float get_bottle_rot_angle() const {return (rotates() ? PI*(0.321*obj_id + color.R + 2.0*color.G) : 0.0);}
	float get_depth () const {return get_length();} // some objects use depth rather than length
	float get_radius() const;
	float get_ortho_radius() const {return (is_hanging() ? get_radius() : 0.5*dz());} // hanging is dim=2/z
	float get_chair_leg_width() const {assert(type == TYPE_CHAIR); return (in_mall() ? CHAIR_LEG_WIDTH_MALL : CHAIR_LEG_WIDTH);}
	cylinder_3dw get_cylinder() const;
	void toggle_lit_state() {flags ^= RO_FLAG_LIT;}
	static bool enable_rugs();
	static bool enable_pictures();
	int get_rug_tid() const;
	int get_picture_tid() const;
	int get_tv_tid() const;
	int get_comp_monitor_tid() const;
	int get_sheet_tid() const;
	int get_paper_tid() const;
	int get_food_box_tid() const;
	std::string const &get_food_box_name() const;
	int get_model_id() const;
	void set_as_bottle(unsigned rand_id, unsigned max_type=NUM_BOTTLE_TYPES-1, bool no_empty=0, unsigned exclude_mask=0);
	void remove() {type = TYPE_BLOCKER; flags = (RO_FLAG_NOCOLL | RO_FLAG_INVIS);} // replace it with an invisible blocker that won't collide with anything
	colorRGBA get_color() const;
	colorRGBA get_model_color() const;
	vector3d get_dir() const {return vector_from_dim_dir(dim, dir);}
	rand_gen_t create_rgen() const;
	ball_type_t const &get_ball_type() const;
}; // room_object_t
typedef vector<room_object_t> vect_room_object_t;

inline void set_obj_id(vect_room_object_t &objs) {objs.back().obj_id = (uint16_t)objs.size();}

struct carried_item_t : public room_object_t {
	unsigned use_count=0;
	carried_item_t() {}
	carried_item_t(room_object_t const &o) : room_object_t(o) {}
	float get_remaining_capacity_ratio() const;
};

struct custom_item_t {
	std::string name;
	float value=0.0, weight=0.0, health=0.0;
	custom_item_t() {}
	custom_item_t(std::string const &n, float v, float w, float h=0.0) : name(n), value(v), weight(w), health(h) {}
	bool valid() const {return !name.empty();}
};
struct room_obj_or_custom_item_t {
	room_object_t obj;
	custom_item_t item;
};

struct room_obj_dstate_t { // state used for dynamic room objects
	vector3d velocity;
	xform_matrix rot_matrix;
};

class vbo_cache_t {
	struct vbo_cache_entry_t {
		unsigned vbo, size;
		vbo_cache_entry_t(unsigned vbo_, unsigned size_) : vbo(vbo_), size(size_) {}
	};
private:
	vector<vbo_cache_entry_t> entries[2]; // {vertex, index}
	unsigned v_alloc=0, v_used=0, v_reuse=0, v_free=0;
	size_t s_alloc=0, s_used=0, s_reuse=0, s_free=0;
public:
	vbo_cache_entry_t alloc( unsigned size, bool is_index=0);
	void free(unsigned &vbo, unsigned size, bool is_index=0);
	void clear(); // unused
	unsigned size() const {return (entries[0].size() + entries[1].size());}
	void print_stats() const;
};

struct rgeom_storage_t {
	typedef vert_norm_comp_tc_color vertex_t;
	typedef vector<vertex_t> vect_vertex_t;
	vector<vertex_t> quad_verts, itri_verts;
	vector<unsigned> indices; // for indexed_quad_verts
	tid_nm_pair_t tex; // used sort of like a map key

	rgeom_storage_t() {}
	rgeom_storage_t(tid_nm_pair_t const &tex_) : tex(tex_) {}
	bool empty() const {return (quad_verts.empty() && itri_verts.empty());}
	void clear(bool free_memory=0);
	void swap_vectors(rgeom_storage_t &s);
	void swap(rgeom_storage_t &s);
	unsigned get_tot_vert_capacity() const {return (quad_verts.capacity() + itri_verts.capacity());}
	unsigned get_mem_usage() const {return (get_cont_mem_usage(quad_verts) + get_cont_mem_usage(itri_verts) + get_cont_mem_usage(indices));}
};

class rgeom_mat_t : public rgeom_storage_t { // simplified version of building_draw_t::draw_block_t

	indexed_vao_manager_with_shadow_t vao_mgr;
	static vbo_cache_t vbo_cache; // shared across all buildings and materials
	cube_t bcube;
public:
	unsigned num_verts=0, num_ixs=0, vert_vbo_sz=0, ixs_vbo_sz=0; // for drawing
	uint8_t dir_mask=0; // {-x, +x, -y, +y, -z, +z}
	bool en_shadows=0;

	rgeom_mat_t(tid_nm_pair_t const &tex_=tid_nm_pair_t()) : rgeom_storage_t(tex_) {}
	//~rgeom_mat_t() {assert(vao_mgr.vbo == 0); assert(vao_mgr.ivbo == 0);} // VBOs should be freed before destruction
	void enable_shadows() {en_shadows = 1;}
	void clear();
	void clear_vbos();
	void clear_vectors(bool free_memory=0) {rgeom_storage_t::clear(free_memory);}
	void add_cube_to_verts(cube_t const &c, colorRGBA const &color, point const &tex_origin=all_zeros, unsigned skip_faces=0,
		bool swap_tex_st=0, bool mirror_x=0, bool mirror_y=0, bool inverted=0, bool z_dim_uses_ty=0, float tx_add=0.0, float ty_add=0.0);
	void add_cube_to_verts_untextured(cube_t const &c, colorRGBA const &color, unsigned skip_faces=0);
	void add_ortho_cylin_to_verts(cube_t const &c, colorRGBA const &color, int dim, bool draw_bot, bool draw_top, bool two_sided=0, bool inv_tb=0,
		float rs_bot=1.0, float rs_top=1.0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0,
		bool swap_txy=0, float len_tc2=1.0, float len_tc1=0.0, int half_or_quarter=0);
	void add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided=0, bool inv_tb=0,
		float rs_bot=1.0, float rs_top=1.0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0,
		bool swap_txy=0, float len_tc2=1.0, float len_tc1=0.0, int half_or_quarter=0);
	void add_vcylin_to_verts_tscale(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top);
	void add_cylin_to_verts(point const &bot, point const &top, float bot_radius, float top_radius, colorRGBA const &color, bool draw_bot, bool draw_top,
		bool two_sided=0, bool inv_tb=0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0,
		bool swap_txy=0, float len_tc2=1.0, float len_tc1=0.0, int half_or_quarter=0);
	void add_disk_to_verts(point const &pos, float radius, vector3d const &dir, colorRGBA const &color, bool swap_txy=0, bool inv_ts=0, bool inv_tt=0, unsigned ndiv=N_CYL_SIDES);
	void add_vert_disk_to_verts(point const &pos, float radius, bool normal_z_neg, colorRGBA const &color, bool swap_txy=0, bool inv_ts=0, bool inv_tt=0, unsigned ndiv=N_CYL_SIDES) {
		add_disk_to_verts(pos, radius, (normal_z_neg ? -plus_z : plus_z), color, swap_txy, inv_ts, inv_tt, ndiv);
	}
	void add_sphere_to_verts(point const &center, vector3d const &size, colorRGBA const &color, bool low_detail=0,
		vector3d const &skip_hemi_dir=zero_vector, tex_range_t const &tr=tex_range_t(), xform_matrix const *const matrix=nullptr, float ts_add=0.0, float tt_add=0.0);
	void add_sphere_to_verts(point const &center, float radius, colorRGBA const &color, bool low_detail=0, tex_range_t const &tr=tex_range_t()) {
		add_sphere_to_verts(center, vector3d(radius, radius, radius), color, low_detail, zero_vector, tr);
	}
	void add_sphere_to_verts(cube_t const &c, colorRGBA const &color, bool low_detail=0, vector3d const &skip_hemi_dir=zero_vector,
		tex_range_t const &tr=tex_range_t(), xform_matrix const *const matrix=nullptr, float ts_add=0.0, float tt_add=0.0)
	{
		add_sphere_to_verts(c.get_cube_center(), 0.5*c.get_size(), color, low_detail, skip_hemi_dir, tr, matrix, ts_add, tt_add);
	}
	void add_vert_torus_to_verts (point const &center, float r_inner, float r_outer, colorRGBA const &color,
		float tscale=1.0, bool low_detail=0, int half_or_quarter=0, float s_offset=0.0, unsigned ndivo=0, unsigned ndivi=0, float spiral_offset=0.0);
	void add_ortho_torus_to_verts(point const &center, float r_inner, float r_outer, unsigned dim, colorRGBA const &color,
		float tscale=1.0, bool low_detail=0, int half_or_quarter=0, float s_offset=0.0, unsigned ndivo=0, unsigned ndivi=0, float spiral_offset=0.0);
	void add_contained_vert_torus_to_verts(cube_t const &c, colorRGBA const &color, float tscale=1.0, bool low_detail=0);
	void add_triangle_to_verts(point const v[3], colorRGBA const &color, bool two_sided, float tscale=1.0);
	void add_quad_to_verts(point const v[4], colorRGBA const &color, float tscale=1.0);
	void create_vbo(building_t const &building);
	void create_vbo_inner();
	void vao_setup(bool shadow_only);
	void draw(tid_nm_pair_dstate_t &state, brg_batch_draw_t *bbd, int shadow_only, bool reflection_pass, bool exterior_geom);
	void pre_draw(int shadow_only) const;
	void draw_geom() const;
	void draw_inner(int shadow_only) const;
	void upload_draw_and_clear(tid_nm_pair_dstate_t &state);
}; // rgeom_mat_t

struct building_materials_t : public vector<rgeom_mat_t> {
	bool valid=0;
	void clear();
	void invalidate() {valid = 0;}
	unsigned count_all_verts() const;
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows);
	void create_vbos(building_t const &building);
	void draw(brg_batch_draw_t *bbd, shader_t &s, int shadow_only, bool reflection_pass, bool exterior_geom=0);
	void upload_draw_and_clear(shader_t &s);
};

struct obj_model_inst_t {
	bool int_vis_only;
	unsigned obj_id;
	vector3d dir;
  obj_model_inst_t(unsigned oid, vector3d const &d, bool ivo) : int_vis_only(ivo), obj_id(oid), dir(d) {}
};
struct obj_model_inst_with_obj_t : obj_model_inst_t {
	obj_model_inst_with_obj_t(obj_model_inst_t const &i, room_object_t const &obj_) : obj_model_inst_t(i), obj(obj_) {}
	room_object_t obj;
};

class brg_batch_draw_t {
	struct mat_entry_t {
		tid_nm_pair_t tex;
		vector<rgeom_mat_t const *> mats;
		mat_entry_t() {}
		mat_entry_t(rgeom_mat_t const &m) : tex(m.tex) {mats.push_back(&m);}
	};
	struct tile_block_t {
		vector<mat_entry_t> to_draw;
		cube_t bcube;
		tile_block_t(cube_t const &bcube_=cube_t()) : bcube(bcube_) {}
	};
	vector<mat_entry_t> to_draw; // interior objects
	vector<tile_block_t> ext_by_tile; // exterior objects, stored by tile for shadow mapping
	vector<int> tid_to_first_mat_map; // -1 is unset
	unsigned cur_tile_slot=0;

	void draw_and_clear_batch(vector<mat_entry_t> &batch, tid_nm_pair_dstate_t &state);
public:
	uint8_t camera_dir_mask=0;
	vector<obj_model_inst_with_obj_t> models_to_draw; // models on building exteriors to draw after buildings

	bool has_ext_geom() const;
	void clear();
	void set_camera_dir_mask(point const &camera_bs, cube_t const &bcube);
	void next_tile(cube_t const &bcube);
	void add_material(rgeom_mat_t const &m, bool is_ext_tile=0);
	void draw_and_clear(shader_t &s);
	void draw_and_clear_ext_tiles(shader_t &s, vector3d const &xlate);
	void clear_ext_tiles();
	void clear_obj_models() {models_to_draw.clear();}
	void draw_obj_models(shader_t &s, vector3d const &xlate, bool shadow_only) const;
};

struct tape_quad_batch_draw : public quad_batch_draw {
	int moving_vert_cyilin_int_tape(point &cur_pos, point const &prev_pos, float z1, float z2, float radius, float slow_amt, bool is_player) const;
	void split_tape_at(unsigned first_vert, point const &pos, float min_zval);
};

struct paint_draw_t {
	quad_batch_draw sp_qbd[NUM_SP_EMISSIVE_COLORS+1], m_qbd; // {spraypaint, markers}
	bool have_any_sp() const;
	quad_batch_draw &get_paint_qbd(bool is_marker, unsigned emissive_color_id);
	void draw_paint(shader_t &s) const;
	void clear();
};
struct building_decal_manager_t {
	paint_draw_t paint_draw[2]; // {interior, exterior}
	quad_batch_draw blood_qbd[2], tp_qbd, pend_tape_qbd, glass_qbd, burn_qbd; // blood_qbd: {red human blood, bug guts or stains}
	tape_quad_batch_draw tape_qbd; // for tape, but not pend_tape because it hasn't been placed yet
	rand_gen_t rgen;

	void commit_pend_tape_qbd();
	void add_burn_spot(point const &pos, float radius);
	void add_blood_or_stain(point const &pos, float radius, colorRGBA const &color, bool is_blood=0, unsigned dim=2, bool dir=1);
	void draw_building_interior_decals(shader_t &s, bool player_in_building, bool shadow_only) const;
};

class particle_manager_t {
	struct particle_t {
		point pos;
		vector3d vel;
		colorRGBA color;
		float init_radius=0.0, radius=0.0, time=0.0, alpha=1.0, coll_radius=1.0;
		int parent_obj_id=0;
		unsigned effect=PART_EFFECT_NONE, bounce_count=0;

		particle_t() {}
		particle_t(point const &p, vector3d const &v, colorRGBA const &c, float r, unsigned e, int pid=-1, float cr=1.0) :
			pos(p), vel(v), color(c), init_radius(r), radius(r), alpha(color.A), coll_radius(cr), parent_obj_id(pid), effect(e) {}
	};
	vector<particle_t> particles;
	quad_batch_draw qbds[NUM_PART_EFFECTS]; // one per perticle effect
	rand_gen_t rgen;
public:
	void add_particle(point const &pos, vector3d const &vel, colorRGBA const &color, float radius, unsigned effect, int pid=-1, float coll_radius=1.0) {
		particles.emplace_back(pos, vel, color, radius, effect, pid, coll_radius);
	}
	void add_for_obj(room_object_t &obj, float pradius, vector3d const &dir, float part_vel, unsigned min_parts, unsigned max_parts, unsigned effect, int parent_obj_id);
	cube_t get_bcube() const;
	void next_frame(building_t &building);
	void add_lights(vector3d const &xlate, building_t const &building, occlusion_checker_noncity_t const &oc, cube_t &lights_bcube) const;
	void draw(shader_t &s, vector3d const &xlate);
};

template<typename T> cube_t get_cube_height_radius(point const &center, T radius, float height) { // T can be float or vector3d
	cube_t c(center);
	c.expand_by_xy(radius);
	c.z2() += height;
	return c;
}

class fire_manager_t {
	struct fire_t {
		point pos; // pos is the bottom
		float max_radius=0.0, radius=0.0, time=0.0, next_smoke_time=0.0;
		fire_t(point const &pos_, float max_radius_) : pos(pos_), max_radius(max_radius_) {}
		float get_height() const {return 4.0*radius;}
		point get_center() const {return pos + vector3d(0.0, 0.0, 0.5*get_height());}
		cube_t get_bcube() const {return get_cube_height_radius(pos, radius, get_height());}
	};
	vector<fire_t> fires;
	quad_batch_draw qbd;
	rand_gen_t rgen;
public:
	void spawn_fire(point const &pos, float size);
	cube_t get_bcube() const;
	bool get_closest_fire(point const &pos, float xy_radius, float z1, float z2, point *fire_pos=nullptr) const;
	void add_fire_bcubes_for_cube(cube_t const &sel_cube, vect_cube_t &fire_bcubes) const;
	void put_out_fires(point const &p1, point const &p2, float radius);
	void next_frame(particle_manager_t &particle_manager);
	void add_lights(vector3d const &xlate, building_t const &building, occlusion_checker_noncity_t const &oc, cube_t &lights_bcube) const;
	void draw(shader_t &s, vector3d const &xlate);
};

struct droplet_spawner_t {
	point pos;
	float radius, period, last_spawned=0.0;
	droplet_spawner_t(point const &pos_, float radius_, float period_) : pos(pos_), radius(radius_), period(period_) {}
};

struct index_pair_t {
	unsigned ix[2] = {};
	index_pair_t() {}
	index_pair_t(vect_cube_t const V[2]) {ix[0] = V[0].size(); ix[1] = V[1].size();}
};

struct courtyard_t : public cube_t {
	int16_t room_ix=-1, door_ix=-1; // starts as <unset>
};

struct door_rotation_t {
	float angle=0.0, shift=0.0;
};
struct door_handle_t {
	point center;
	float height;
	bool dim, hdir, mirror, closed;
	vector3d dir;
	door_handle_t(point const &cp, float h, bool D, bool hd, bool m, bool c, vector3d const d) : center(cp), height(h), dim(D), hdir(hd), mirror(m), closed(c), dir(d) {}
};

struct building_room_geom_t {

	bool has_pictures=0, has_garage_car=0, modified_by_player=0, have_clock=0, glass_floor_split=0;
	unsigned char num_pic_tids=0, invalidate_mats_mask=0;
	float obj_scale=1.0;
	unsigned wall_ps_start=0, buttons_start=0, stairs_start=0, backrooms_start=0, retail_start=0; // index of first object of {TYPE_PG_*|TYPE_PSPACE, TYPE_BUTTON, TYPE_STAIR, retail}
	unsigned init_num_doors=0, init_num_dstacks=0; // required for removing doors added by backrooms generation when room_geom is deleted
	unsigned pool_ramp_obj_ix=0, pool_stairs_start_ix=0, last_animal_update_frame=0;
	point tex_origin;
	colorRGBA wood_color;
	courtyard_t courtyard;
	// objects in rooms; expanded_objs is for things that have been expanded for player interaction; model_objs is for models in drawers; trim_objs is for wall/door/window trim
	vect_room_object_t objs, expanded_objs, model_objs, trim_objs;
	vector<room_obj_dstate_t> obj_dstate;
	vector<obj_model_inst_t> obj_model_insts;
	vector<door_handle_t> door_handles; // for 3D model drawing
	vector<unsigned> moved_obj_ids;
	vect_rat_t    rats, sewer_rats, pet_rats;
	vect_spider_t spiders, sewer_spiders;
	vect_snake_t  snakes, pet_snakes;
	vect_bird_t   pet_birds;
	vect_insect_t insects;
	// {large static, small static, dynamic, lights, alpha mask, transparent, door} materials
	building_materials_t mats_static, mats_small, mats_text, mats_detail, mats_dynamic, mats_lights, mats_amask, mats_alpha, mats_doors, mats_exterior, mats_ext_detail;
	rgeom_mat_t mats_glass[2]; // {viewed from below, viewed from above}
	vect_cube_t light_bcubes, shelf_rack_occluders[2], glass_floors; // shelf_rack_occluders: {back, top}
	vect_cube_t pgbr_walls[2]; // parking garage and backrooms walls, in each dim
	vector<index_pair_t> pgbr_wall_ixs; // indexes into pgbr_walls
	building_decal_manager_t decal_manager;
	particle_manager_t particle_manager;
	fire_manager_t fire_manager;
	vector<droplet_spawner_t> droplet_spanwers; // for flooded basements

	building_room_geom_t(point const &tex_origin_=all_zeros) : tex_origin(tex_origin_), wood_color(WHITE) {}
	bool empty() const {return objs.empty();}
	void clear();
	void clear_materials();
	void invalidate_static_geom  () {invalidate_mats_mask |= (1 << MAT_TYPE_STATIC );}
	void invalidate_model_geom   () {invalidate_static_geom();}
	void invalidate_small_geom   () {invalidate_mats_mask |= (1 << MAT_TYPE_SMALL  );}
	void update_text_draw_data   () {invalidate_mats_mask |= (1 << MAT_TYPE_TEXT   );}
	void invalidate_lights_geom  () {invalidate_mats_mask |= (1 << MAT_TYPE_LIGHTS );} // cache state and apply change later in case this is called from a different thread
	void invalidate_detail_geom  () {invalidate_mats_mask |= (1 << MAT_TYPE_DETAIL );}
	void update_dynamic_draw_data() {invalidate_mats_mask |= (1 << MAT_TYPE_DYNAMIC);}
	void invalidate_door_geom    () {invalidate_mats_mask |= (1 << MAT_TYPE_DOORS  );}
	void check_invalid_draw_data();
	void invalidate_draw_data_for_obj(room_object_t const &obj, bool was_taken=0);
	unsigned get_num_verts() const;
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool transparent=0, bool exterior=0) {
		return get_building_mat(tex, dynamic, small, transparent, exterior).get_material(tex, inc_shadows);
	}
	rgeom_mat_t &get_untextured_material(bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool transparent=0, bool exterior=0) {
		return get_material(tid_nm_pair_t(-1, 1.0, inc_shadows, transparent), inc_shadows, dynamic, small, transparent, exterior);
	}
	rgeom_mat_t &get_wood_material(float tscale=1.0, bool inc_shadows=1, bool dynamic=0, unsigned small=0, bool exterior=0);
	rgeom_mat_t &get_metal_material(bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool exterior=0, colorRGBA const &spec_color=WHITE);
	rgeom_mat_t &get_scratched_metal_material(float tscale, bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool exterior=0);
	colorRGBA apply_wood_light_color(room_object_t const &o) const;
	void add_tquad(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
		colorRGBA const &color, bool invert_tc_x, bool exclude_frame, bool no_tc);
	vect_room_object_t::const_iterator get_placed_objs_end() const {return (objs.begin() + buttons_start);} // excludes buttons, stairs, and elevators
	vect_room_object_t::const_iterator get_stairs_start   () const {return (objs.begin() + stairs_start );} // excludes stairs
	bool cube_int_backrooms_walls(cube_t const &c) const;
	// Note: these functions are all for drawing objects / adding them to the vertex list
	void add_tc_legs(cube_t const &c, room_object_t const &obj, colorRGBA const &color, float width, bool recessed, float tscale,
		bool use_metal_mat=0, bool draw_tops=0, float frame_height=0.0);
	void add_table(room_object_t const &c, float tscale, float top_dz, float leg_width);
	void add_chair(room_object_t const &c, float tscale);
	void add_dresser(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void draw_mirror_surface(room_object_t const &c, cube_t const &mirror, bool dim, bool dir, bool shadowed);
	void add_dresser_mirror(room_object_t const &c, float tscale);
	void add_dresser_drawers(room_object_t const &c, float tscale);
	void add_drawers(room_object_t const &c, float tscale, vect_cube_t const &drawers, unsigned drawer_index_offset=0);
	void add_stair(room_object_t const &c, float tscale, vector3d const &tex_origin, bool is_small_pass);
	void add_stairs_wall(room_object_t const &c, vector3d const &tex_origin, tid_nm_pair_t const &wall_tex);
	void add_wall_or_pillar (room_object_t const &c, vector3d const &tex_origin, tid_nm_pair_t const &wall_tex);
	void add_basement_pillar(room_object_t const &c, tid_nm_pair_t const &wall_tex);
	void add_basement_beam  (room_object_t const &c, tid_nm_pair_t const &wall_tex);
	void add_shelf_wall     (room_object_t const &c, tid_nm_pair_t const &wall_tex);
	void add_parking_space(room_object_t const &c, float tscale);
	void add_ramp(room_object_t const &c, float thickness, bool skip_bottom, rgeom_mat_t &mat);
	void add_pg_ramp(room_object_t const &c, float tscale);
	void add_pipe(room_object_t const &c, bool add_exterior);
	void add_sprinkler(room_object_t const &c);
	void add_valve(room_object_t const &c);
	void add_gauge(room_object_t const &c);
	void add_duct(room_object_t const &c);
	void add_warning_light(room_object_t const &c);
	void add_pallet(room_object_t const &c);
	void add_curb(room_object_t const &c);
	void add_chimney(room_object_t const &c, tid_nm_pair_t const &tex);
	void add_breaker_panel(room_object_t const &c);
	void add_attic_door(room_object_t const &c, float tscale);
	void add_attic_interior_and_rafters(building_t const &b, float tscale, bool detail_pass);
	void add_skylights_details(building_t const &b);
	void add_skylight_details(cube_t const &skylight, bool has_skylight_light);
	void add_elevator(room_object_t const &c, elevator_t const &e, float tscale, float fc_thick_scale,
		unsigned floor_offset, float floor_spacing, float window_vspace, bool has_parking_garage, bool is_powered);
	void add_escalator(escalator_t const &e, float floor_spacing, bool draw_static, bool draw_dynamic);
	void add_tunnel(tunnel_seg_t const &t);
	void add_tunnel_water(tunnel_seg_t const &t);
	void add_elevator_doors(elevator_t const &e, float fc_thick_scale);
	void add_light(room_object_t const &c, float tscale);
	void add_rug(room_object_t const &c);
	void add_blanket(room_object_t const &c);
	void add_picture(room_object_t const &c);
	static void add_book_title(std::string const &title, cube_t const &title_area, rgeom_mat_t &mat, colorRGBA const &color,
		unsigned hdim, unsigned tdim, unsigned wdim, bool cdir, bool ldir, bool wdir);
	void add_book(room_object_t const &c, bool inc_lg, bool inc_sm, bool inc_text, float tilt_angle=0.0, unsigned extra_skip_faces=0, bool no_title=0, float z_rot_angle=0.0);
	void add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, bool inc_text, float tscale,
		bool no_shelves=0, float sides_scale=1.0, point const *const use_this_tex_origin=nullptr, vect_room_object_t *books=nullptr);
	void add_wine_rack(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale);
	void add_desk(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void add_reception_desk(room_object_t const &c, float tscale);
	void add_conference_table(room_object_t const &c, float tscale);
	void add_bed(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale);
	void add_window(room_object_t const &c, float tscale);
	void add_int_window(room_object_t const &c);
	void add_crack(room_object_t const &c);
	void add_switch(room_object_t const &c, bool draw_detail_pass);
	void add_breaker(room_object_t const &c);
	void add_flat_textured_detail_wall_object(room_object_t const &c, colorRGBA const &side_color, int tid, bool skip_z1_face,
		bool draw_all_faces=0, bool detail=1, bool mirror_y=0);
	void add_outlet(room_object_t const &c);
	void add_vent(room_object_t const &c);
	void add_plate(room_object_t const &c);
	void add_water_plane(room_object_t const &c, cube_t const &water_area, float water_level);
	void add_tub_outer (room_object_t const &c);
	void add_sink_water(room_object_t const &c);
	void add_tv_picture(room_object_t const &c);
	void add_cup_liquid(room_object_t const &c);
	void add_trashcan(room_object_t const &c);
	void add_bucket(room_object_t const &c, bool draw_metal, bool draw_liquid);
	void add_laundry_basket(room_object_t const &c);
	void add_water_heater(room_object_t const &c);
	void add_furnace(room_object_t const &c);
	void add_server(room_object_t const &c);
	void add_pool_ball(room_object_t const &c);
	void add_pool_cue (room_object_t const &c);
	void add_wall_mount(room_object_t const &c);
	void add_toaster_proxy(room_object_t const &c);
	void add_br_stall(room_object_t const &c);
	void add_cubicle(room_object_t const &c, float tscale);
	void add_window_sill(room_object_t const &c);
	void add_exterior_step(room_object_t const &c);
	void add_balcony(room_object_t const &c, float ground_floor_z1, bool is_in_city);
	void add_sign(room_object_t const &c, bool inc_back, bool inc_text, bool exterior_only=0);
	void add_false_door_int(room_object_t const &c);
	void add_false_door(room_object_t const &c);
	void add_counter(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void add_cabinet(room_object_t const &c, room_object_t const &parent, float tscale, bool inc_lg, bool inc_sm);
	void add_dishwasher_front(room_object_t const &c, float back_wall, unsigned dw_skip_faces);
	void add_dishwasher(room_object_t const &c);
	void add_closet(room_object_t const &c, tid_nm_pair_t const &wall_tex, colorRGBA const &trim_color, bool inc_lg, bool inc_sm);
	void add_hanger_rod(room_object_t const &c);
	void add_drain_cover(cube_t const &c, colorRGBA const &color);
	void add_drain_pipe(room_object_t const &c);
	void add_electrical_wire(room_object_t const &c, vector3d const &rot_axis);
	void add_electrical_wire_pair(room_object_t const &c);
	void add_key(room_object_t const &c);
	void add_money(room_object_t const &c);
	void add_phone(room_object_t const &c);
	void add_tproll(room_object_t const &c);
	void add_tape(room_object_t const &c);
	static void add_spraycan_to_material(room_object_t const &c, rgeom_mat_t &mat, bool draw_bottom=0);
	void add_spraycan(room_object_t const &c);
	void add_button(room_object_t const &c, bool inc_geom, bool inc_text);
	void add_crate(room_object_t const &c);
	void add_box(room_object_t const &c);
	void add_paint_can(room_object_t const &c);
	void add_shelves(room_object_t const &c, float tscale);
	void add_rack(room_object_t const &c, bool add_rack, bool add_objs, bool obj_text_pass=0);
	void add_chimney_cap(room_object_t const &c);
	void add_ext_ladder(room_object_t const &c);
	void add_int_ladder(room_object_t const &c);
	void add_catwalk(room_object_t const &c);
	void add_machine_pipe_in_region(room_object_t const &c, cube_t const &region, float rmax, unsigned dim,
		vect_sphere_t &pipe_ends, rand_gen_t &rgen, bool allow_valve=0, bool valve_on_cylin=0, bool add_coil=0, bool is_high=0);
	void add_spring(point pos, float radius, float r_wire, float length, float coil_gap, unsigned dim,
		colorRGBA const &color, colorRGBA const &spec_color=WHITE, bool sparse=0);
	void add_machine(room_object_t const &c, float floor_ceil_gap, bldg_industrial_info_t const *ind_info);
	void add_keyboard(room_object_t const &c);
	void add_obj_with_top_texture  (room_object_t const &c, std::string const &texture_name, colorRGBA const &sides_color, bool is_small=0);
	void add_obj_with_front_texture(room_object_t const &c, std::string const &texture_name, colorRGBA const &front_color, colorRGBA const &sides_color, bool is_small=0);
	void add_obj_with_front_texture(room_object_t const &c, std::string const &texture_name, colorRGBA const &sides_color, bool is_small=0) {
		add_obj_with_front_texture(c, texture_name, c.color, sides_color, is_small); // front_color = c.color
	}
	void add_computer(room_object_t const &c);
	void add_laptop(room_object_t const &c);
	void add_pizza_box(room_object_t const &c);
	void add_pizza_top(room_object_t const &c);
	void add_mwave(room_object_t const &c);
	void add_vending_machine(room_object_t const &c);
	void add_cabinet_with_open_door(room_object_t const &c, cube_t const &door, colorRGBA const &color,
		float wall_thickness, unsigned front_face_mask, float z1_adj=0.0, float back_adj=0.0);
	void add_locker(room_object_t const &c);
	void add_mirror (room_object_t const &c);
	void add_med_cab(room_object_t const &c);
	rgeom_mat_t &get_shower_tile_mat(room_object_t const &c, float tscale, colorRGBA &color);
	void draw_shower_head(room_object_t const &shower, float radius, float wall_pos, float extent_amt, bool head_dim);
	void add_shower(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void add_shower_tub(room_object_t const &c, tid_nm_pair_t const &wall_tex, colorRGBA const &trim_color, float tscale, bool inc_lg, bool inc_sm);
	void add_bottle(room_object_t const &c, bool add_bottom=0, float label_rot_angle=0.0);
	void add_drink_can(room_object_t const &c, bool add_bottom=0);
	void add_vase(room_object_t const &c);
	void add_paper(room_object_t const &c);
	static void add_pen_pencil_marker_to_material(room_object_t const &c_, rgeom_mat_t &mat);
	void add_pen_pencil_marker(room_object_t const &c);
	void add_flooring (room_object_t const &c, float tscale);
	void add_pool_tile(room_object_t const &c, float tscale);
	void add_wall_trim(room_object_t const &c, bool for_closet=0);
	void add_blinds(room_object_t const &c);
	void add_fireplace(room_object_t const &c, float tscale);
	void add_filing_cabinet(room_object_t const &c, bool inc_lg, bool inc_sm);
	void add_stapler(room_object_t const &c);
	void add_eraser(room_object_t const &c);
	void add_fire_ext_mount(room_object_t const &c);
	void add_fire_ext_sign (room_object_t const &c);
	void add_teeshirt(room_object_t const &c);
	void add_pants(room_object_t const &c);
	void add_ceiling_fan_light(room_object_t const &fan, room_object_t const &light);
	void add_railing(room_object_t const &c);
	void add_downspout(room_object_t const &c);
	void add_plant_pot(room_object_t const &c, float height, float radius, bool no_dirt);
	void add_potted_plant(room_object_t const &c, bool inc_pot, bool inc_plant);
	void add_tree(room_object_t const &c, bool inc_pot, bool inc_tree);
	void add_lg_ball(room_object_t const &c);
	void add_toy(room_object_t const &c);
	void add_pan(room_object_t const &c);
	void add_pool_float(room_object_t const &c);
	void add_bench(room_object_t const &c);
	void add_diving_board(room_object_t const &c);
	void add_flashlight_to_material(room_object_t const &c, rgeom_mat_t &mat, unsigned ndiv);
	void add_flashlight(room_object_t const &c);
	void add_candle(room_object_t const &c);
	void add_camera(room_object_t const &c);
	void add_clock (room_object_t const &c, bool add_dynamic);
	void add_food_box(room_object_t const &c);
	void add_safe(room_object_t const &c);
	void add_checkout(room_object_t const &c, float tscale);
	void add_fishtank(room_object_t const &c);
	void add_metal_bar(room_object_t const &c);
	void add_ibeam(room_object_t const &c);
	void add_chem_tank(room_object_t const &c, bool draw_label);
	void add_hvac_unit(room_object_t const &c);
	void add_vent_fan_frame(room_object_t const &c);
	void add_store_gate(cube_t const &c, bool dim, float open_amt);
	void add_theft_sensor(room_object_t const &c, bool alarm_mode=0);
	void add_lava_lamp(room_object_t const &c);
	void add_trash(room_object_t const &c);
	void add_spider_web(room_object_t const &c);
	void add_pet_cage(room_object_t const &c);
	void add_debug_shape(room_object_t const &c);
	static void draw_ball_in_building(room_object_t  const &c, shader_t &s);
	void draw_interactive_player_obj(carried_item_t const &c, shader_t &s, vector3d const &xlate);
	// functions for expanding nested objects
	void expand_shelves(room_object_t const &c, bool add_models_mode=0);
	void expand_shelfrack(room_object_t const &c);
	void get_bookcase_books(room_object_t const &c, vect_room_object_t &books) {add_bookcase(c, 0, 0, 0, 1.0, 0, 1.0, nullptr, &books);} // Note: technically const
	void expand_closet(room_object_t const &c) {add_closet_objects(c, expanded_objs);}
	void expand_cabinet(room_object_t const &c);
	void expand_wine_rack(room_object_t const &c) {add_wine_rack_bottles(c, expanded_objs);}
	void expand_med_cab(room_object_t const &c);
	void expand_locker (room_object_t const &c);
	void expand_breaker_panel(room_object_t const &c, building_t const &building);
	void expand_dishwasher(room_object_t &c, cube_t const &dishwasher);
	void unexpand_dishwasher(room_object_t &c, cube_t const &dishwasher);
	bool expand_object(room_object_t &c, building_t const &building);
	static room_object_t get_item_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix, unsigned item_ix, float &stack_z1);
	static void add_draw_items(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix, vect_room_object_t &objects);
	bool maybe_spawn_spider_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_id, float floor_spacing, bool is_door);
	// other functions
	bool closet_light_is_on(cube_t const &closet) const;
	bool is_inside_closed_locker(room_object_t const &obj) const;
	int find_nearest_pickup_object(building_t const &building, point const &at_pos, vector3d const &in_dir, float range, float &obj_dist) const;
	bool cube_intersects_moved_obj(cube_t const &c, int ignore_obj_id=-1) const;
	bool open_nearest_drawer(building_t &building, point const &at_pos, vector3d const &in_dir, float range, bool pickup_item, bool check_only);
	void replace_with_hanging_wires(room_object_t &obj, room_object_t const &old_obj, float wire_radius, bool vertical);
	void remove_object(unsigned obj_id, building_t &building);
	bool player_pickup_object(building_t &building, point const &at_pos, vector3d const &in_dir);
	void set_factory_machine_seed(unsigned rseed1, unsigned rseed2);
	void set_factory_machine_seed_from_obj(room_object_t const &obj) {set_factory_machine_seed(obj.item_flags, obj.obj_id);}
	void update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken);
	room_object_t &get_room_object_by_index(unsigned obj_id);
	float get_combined_obj_weight(room_object_t const &obj) const;
	int find_avail_obj_slot() const;
	void add_expanded_object(room_object_t const &obj);
	bool add_room_object(room_object_t const &obj, building_t &building, bool set_obj_id=0, vector3d const &velocity=zero_vector);
	void draw(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, building_t const &building, occlusion_checker_noncity_t &oc,
		vector3d const &xlate, unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building);
	void draw_animals(shader_t &s, building_t const &building, occlusion_checker_noncity_t const &oc, vector3d const &xlate,
		point const &camera_bs, bool shadow_only, bool reflection_pass, bool check_clip_cube) const;
	unsigned allocate_dynamic_state();
	room_obj_dstate_t &get_dstate(room_object_t const &obj);
private:
	building_materials_t &get_building_mat(tid_nm_pair_t const &tex, bool dynamic, unsigned small, bool transparent, bool exterior);
	void create_static_vbos(building_t const &building);
	void create_small_static_vbos(building_t const &building);
	void create_text_vbos();
	void add_text_objs_to_verts(vect_room_object_t const &objs_to_add);
	void create_detail_vbos(building_t const &building);
	void add_nested_objs_to_verts(vect_room_object_t const &objs_to_add);
	void add_small_static_objs_to_verts(vect_room_object_t const &objs_to_add, colorRGBA const &trim_color,
		bool inc_text=0, float floor_ceil_gap=0.0, bldg_industrial_info_t const *ind_info=nullptr);
	void create_obj_model_insts(building_t const &building);
	void create_lights_vbos(building_t const &building);
	void create_dynamic_vbos(building_t const &building, point const &camera_bs, vector3d const &xlate, bool play_clock_tick);
	void create_door_vbos(building_t const &building);
	void add_door_handle(door_t const &door, door_rotation_t const &drot, colorRGBA const &color, bool residential);
	void maybe_add_door_sign(door_t const &door, door_rotation_t const &drot);
	static void add_closet_objects(room_object_t const &c, vect_room_object_t &objects);
	static void get_shelf_objects(room_object_t const &c_in, cube_t const shelves[MAX_SHELVES], unsigned num_shelves, vect_room_object_t &objects, bool add_models_mode=0);
public:
	static void get_shelfrack_objects(room_object_t const &c, vect_room_object_t &objects, bool add_models_mode=0, bool books_only=0);
	static void add_hangers_and_clothing(float window_vspacing, unsigned num_hangers, unsigned flags, int hanger_model_id, int clothing_model_id,
		vect_room_object_t &objects, rand_gen_t &rgen);
	static void add_pet_container(cube_t const &tank, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool in_pet_store,
		vect_room_object_t &objects, rand_gen_t &rgen, unsigned animal_type=TYPE_FISH, unsigned shelf_ix=0);
private:
	static void add_wine_rack_bottles(room_object_t const &c, vect_room_object_t &objects);
	static void add_vert_roll_to_material(room_object_t const &c, rgeom_mat_t &mat, float sz_ratio=1.0, bool player_held=0);
	void add_bcase_book(room_object_t const &c, cube_t const &book, bool inc_lg, bool inc_sm, bool inc_text, bool backwards, bool in_set,
		unsigned skip_faces, unsigned book_ix, unsigned set_start_ix, colorRGBA const &color, float tilt_angle, vect_room_object_t *books);
	void remove_objs_contained_in(cube_t const &c, vect_room_object_t &obj_vect, building_t &building);
}; // building_room_geom_t

struct elevator_t : public oriented_cube_t { // dim/dir applies to the door
	struct call_request_t {
		unsigned floor_ix;
		float zval;
		uint8_t req_dirs; // 1 bit is down, 2 bit is up
		bool inside_press;
		call_request_t(unsigned f, float z, unsigned d, bool ip) : floor_ix(f), zval(z), req_dirs(d), inside_press(ip) {}
		bool operator<(call_request_t const &cr) const {return (cr.inside_press < inside_press);} // sort so that CRs with inside_press=1 are first
	};
	// at_edge is used to avoid expanding the roof elevator cap
	bool at_edge=0, going_up=0, at_dest=0, stop_on_passing_floor=0, hold_doors=0, hold_movement=0, under_skylight=0;
	bool is_moving=0, interior_room=0, in_backrooms=0, no_buttons=0;
	uint8_t in_mall=0; // 0=not in mall, 1=in mall concourse, 2=in mall back area/hallway
	uint8_t adj_pair_ix=0; // 0=not a pair, 1=first, 2=second
	unsigned room_id=0, car_obj_id=0, light_obj_id=0, button_id_start=0, button_id_end=0, num_occupants=0;
	uint64_t skip_floors_mask=0; // good for up to 64 floors
	int at_dest_frame=0, adj_elevator_ix=-1;
	float open_amt=0;
	deque<call_request_t> call_requests; // used as a queue
	vector<unsigned> all_room_ids;

	elevator_t(cube_t const &c, unsigned rid, bool dim_, bool dir_, bool at_edge_, bool interior_room_, uint8_t in_mall_=0, bool in_br=0) :
		oriented_cube_t(c, dim_, dir_), at_edge(at_edge_), interior_room(interior_room_), in_backrooms(in_br), in_mall(in_mall_), room_id(rid)
	{assert(is_strictly_normalized());}
	float get_wall_thickness () const {return 0.02*get_width();}
	float get_frame_width    () const {return ELEVATOR_FRAME_WIDTH*get_width();} // to each side of elevator door
	unsigned get_door_face_id() const {return (2*dim + dir);}
	bool was_called          () const {return !call_requests.empty();}
	bool may_be_moving       () const {return (is_moving || was_called());} // Note: conservative; includes elevator or doors moving
	unsigned get_target_floor() const {assert(was_called()); return call_requests.front().floor_ix;}
	float    get_target_zval () const {assert(was_called()); return call_requests.front().zval;}
	bool was_floor_called(unsigned floor_ix, unsigned up_down_mask) const;
	bool skip_floor_ix   (unsigned floor_ix) const {return (floor_ix < 64 && (skip_floors_mask & (1ULL << floor_ix)));}
	void set_skip_floor  (unsigned floor_ix) {if (floor_ix < 64) {skip_floors_mask |= (1ULL << floor_ix);}}
	unsigned get_coll_cubes(cube_t cubes[5]) const; // returns 1 or 5 cubes
	void call_elevator(unsigned floor_ix, float targ_z, unsigned req_dirs, bool inside_press);
	void register_at_dest();
	void move_closest_in_dir_to_front(float zval, bool dir);
	cube_t get_bcube_padded(float front_pad, float all_sides_pad=0.0) const;
};

unsigned const NUM_RTYPE_SLOTS = 8; // enough for houses; hard max is 8 so that a uint8_t bit mask can be used
inline unsigned wrap_room_floor(unsigned floor) {return min(floor, NUM_RTYPE_SLOTS-1U);}

// room flags
unsigned const ROOM_FLAG_CSTAIRS  = 0x0001; // has center stairs
unsigned const ROOM_FLAG_OFF_FP   = 0x0002; // has office building floorplan
unsigned const ROOM_FLAG_SKYLIGHT = 0x0004; // has a ceiling skylight
unsigned const ROOM_FLAG_IS_ENTRY = 0x0008; // is unit entryway room
unsigned const ROOM_FLAG_MIRROR   = 0x0010; // contains a mirror
unsigned const ROOM_FLAG_HAS_OOO  = 0x0020; // has out-of-order sign
unsigned const ROOM_FLAG_INT_WIND = 0x0040; // has interior window
unsigned const ROOM_FLAG_RO_GEOM  = 0x0080; // room should not have objects placed in it
unsigned const ROOM_FLAG_TUNNEL   = 0x0100; // room has a tunnel connection at one end
unsigned const ROOM_FLAG_NESTED   = 0x0200; // room is nested inside another room
unsigned const ROOM_FLAG_HAS_SUB  = 0x0400; // room has a sub-room nested inside it

struct room_t : public cube_t { // size=56
	bool is_hallway=0, is_office=0, is_sec_bldg=0, is_single_floor=0;
	uint16_t flags=0;
	uint8_t has_stairs=0; // per-floor bit mask; always set to 255 for stairs that span the entire room
	uint8_t has_elevator=0; // number of elevators, usually either 0 or 1
	uint8_t interior=0; // 0=not interior (has windows), 1=interior, 2=extended basement, {3,4}=extended basement connector, dim=interior-3
	uint8_t ext_sides=0; // sides that have exteriors, and likely windows (bits for x1, x2, y1, y2)
	uint8_t part_id=0, num_lights=0, rtype_locked=0;
	uint8_t unit_id=0; // for apartments and hotels
	uint8_t open_wall_mask=0; // {dim x dir}
	uint8_t floors_pre_assigned=0; // one per NUM_RTYPE_SLOTS
	room_type rtype[NUM_RTYPE_SLOTS]; // this applies to the first few floors because some rooms can have variable per-floor assignment
	uint32_t lit_by_floor=0; // used for AI placement; 32 floor is enough for most buildings
	float light_intensity=0.0; // due to room lights, if turned on

	room_t() {assign_all_to(RTYPE_NOTSET, 0);} // locked=0
	room_t(cube_t const &c, unsigned p, unsigned nl=0, bool is_hallway_=0, bool is_office_=0, bool is_sec_bldg_=0, uint8_t interior_=0);
	room_t(room_t const &r, cube_t const &c) {*this = r; copy_from(c);} // sub-room
	void assign_all_to(room_type rt, bool locked=1); // locked by default
	void assign_to(room_type rt, unsigned floor, bool locked=0); // unlocked by default
	void mark_open_wall     (bool dim, bool dir) {open_wall_mask |= (1 << (2*dim + dir));}
	void mark_open_wall_dim (bool dim          ) {open_wall_mask |= (1 << 2*dim) | (1 << (2*dim + 1));} // both dirs for this dim
	bool has_open_wall      (bool dim, bool dir) const {return (open_wall_mask & (1 << (2*dim + dir)));}
	room_type get_room_type (unsigned floor) const {return rtype[wrap_room_floor(floor)];}
	bool is_rtype_locked    (unsigned floor) const {return (rtype_locked & (1 << wrap_room_floor(floor)));}
	bool is_lit_on_floor    (unsigned floor) const {return (lit_by_floor & (1ULL << (floor&31)));}
	bool has_stairs_on_floor(unsigned floor) const {return (get_has_center_stairs() || (has_stairs & (1U << min(floor, 7U))));} // floors >= 7 are treated as the top floor
	bool is_garage_or_shed  (unsigned floor) const {return (is_sec_bldg || get_room_type(floor) == RTYPE_GARAGE || get_room_type(floor) == RTYPE_SHED);}
	bool is_ext_basement     () const {return (interior >= 2);}
	bool is_ext_basement_conn() const {return (interior >= 3);}
	bool inc_half_walls      () const {return (is_hallway || get_office_floorplan() || is_ext_basement());} // hallway, office, or extended basement
	bool is_parking          () const {return (get_room_type(0) == RTYPE_PARKING  );}
	bool is_backrooms        () const {return (get_room_type(0) == RTYPE_BACKROOMS);}
	bool is_mall             () const {return (get_room_type(0) == RTYPE_MALL     );}
	bool is_store            () const {return (get_room_type(0) == RTYPE_STORE    );}
	bool is_retail           () const {return (get_room_type(0) == RTYPE_RETAIL   );}
	bool is_factory          () const {return (get_room_type(0) == RTYPE_FACTORY  );}
	bool is_warehouse        () const {return (get_room_type(0) == RTYPE_WAREHOUSE);}
	bool is_industrial       () const {return (is_factory() || is_warehouse());}
	bool is_mall_or_store    () const {return (is_mall() || is_store());}
	bool is_single_large_room() const {return(is_parking() || is_backrooms() || is_retail() || is_mall() || is_industrial());}
	bool is_single_large_room_or_store() const {return (is_single_large_room() || is_store());}
	bool is_apt_or_hotel_room() const {return (unit_id > 0);}
	bool has_room_of_type(room_type type) const;
	void init_pre_populate(bool is_first_pass);
	float get_light_amt() const;
	unsigned get_floor_containing_zval(float zval, float floor_spacing) const {return (is_single_floor ? 0 : unsigned((zval - z1())/floor_spacing));}
	room_type get_room_type_for_zval  (float zval, float floor_spacing) const {return get_room_type(get_floor_containing_zval(zval, floor_spacing));}

	void set_has_center_stairs() {flags |= ROOM_FLAG_CSTAIRS ;}
	void set_office_floorplan () {flags |= ROOM_FLAG_OFF_FP  ;}
	void set_has_skylight     () {flags |= ROOM_FLAG_SKYLIGHT;}
	void set_is_entryway      () {flags |= ROOM_FLAG_IS_ENTRY;}
	void set_has_mirror       () {flags |= ROOM_FLAG_MIRROR  ;}
	void set_has_out_of_order () {flags |= ROOM_FLAG_HAS_OOO ;}
	void set_interior_window  () {flags |= ROOM_FLAG_INT_WIND;}
	void set_no_geom          () {flags |= ROOM_FLAG_RO_GEOM ;}
	void set_has_tunnel_conn  () {flags |= ROOM_FLAG_TUNNEL  ;}
	void set_is_nested        () {flags |= ROOM_FLAG_NESTED  ;}
	void set_has_subroom      () {flags |= ROOM_FLAG_HAS_SUB ;}
	bool get_has_center_stairs() const {return (flags & ROOM_FLAG_CSTAIRS );}
	bool get_office_floorplan () const {return (flags & ROOM_FLAG_OFF_FP  );}
	bool get_has_skylight     () const {return (flags & ROOM_FLAG_SKYLIGHT);}
	bool get_is_entryway      () const {return (flags & ROOM_FLAG_IS_ENTRY);}
	bool get_has_mirror       () const {return (flags & ROOM_FLAG_MIRROR  );}
	bool get_has_out_of_order () const {return (flags & ROOM_FLAG_HAS_OOO );}
	bool has_interior_window  () const {return (flags & ROOM_FLAG_INT_WIND);}
	bool has_no_geom          () const {return (flags & ROOM_FLAG_RO_GEOM );}
	bool has_tunnel_conn      () const {return (flags & ROOM_FLAG_TUNNEL  );}
	bool is_nested            () const {return (flags & ROOM_FLAG_NESTED  );}
	bool has_subroom          () const {return (flags & ROOM_FLAG_HAS_SUB );}
}; // room_t

struct extb_room_t : public cube_t { // extended basement room candidate
	cube_t conn_bcube;
	bool is_hallway, has_stairs, hallway_dim, connect_dir; // Note: hallway_dim and connect_dir only used for connecting between two buildings
	extb_room_t(cube_t const &c, bool is_hallway_=0, bool has_stairs_=0, bool dim=0, bool dir=0) :
		cube_t(c), is_hallway(is_hallway_), has_stairs(has_stairs_), hallway_dim(dim), connect_dir(dir) {}
	void clip_hallway_to_conn_bcube(bool dim);
};
typedef vector<extb_room_t> vect_extb_room_t;

struct store_info_t {
	bool dim, dir;
	unsigned room_id, store_type, item_category;
	colorRGBA logo_color;
	std::string name;
	vect_cube_t occluders;

	store_info_t(bool d, bool D, unsigned rid, unsigned type, unsigned cat, colorRGBA const &color, std::string const &n) :
		dim(d), dir(D), room_id(rid), store_type(type), item_category(cat), logo_color(color), name(n) {}
	std::string get_full_name() const;
};

struct tunnel_conn_t {
	unsigned dim;
	bool dir;
	float pos, length, radius, water_level=0.0, water_flow=0.0;
	tunnel_conn_t(unsigned dim_, bool dir_, float p, float len, float r) : dim(dim_), dir(dir_), pos(p), length(len), radius(r) {}
};
struct tunnel_seg_t {
	bool dim=0, room_conn=0, room_dir=0, has_gate=0, conns_added=0, closed_ends[2]={};
	int conn_ix[2]={-1,-1}; // index of adjacent connected tunnel in each dir; -1 is none
	unsigned conn_room_ix=0, tseg_ix=0;
	int8_t add_bend_dir[2] = { -1, -1 }; // one per end; -1 is no bend
	float radius=0.0, gate_pos=0.0, water_level=0.0, water_flow=0.0;
	point p[2];
	cube_t bcube, bcube_ext, bcube_draw; // bcube_ext includes the area connecting to the door when room_conn=1 and overlap area at bends
	vector<tunnel_conn_t> conns; // connections to other tunnels

	tunnel_seg_t(point const &p1, point const &p2, float radius_);
	float get_length    () const {return (p[1][dim] - p[0][dim]);}
	float get_bar_radius() const {return 0.025*radius;}
	void set_gate(float pos) {gate_pos = pos; has_gate = 1;}
	void set_as_room_conn(bool rdir, float wall_gap);
	cube_t get_walk_area(point const &pos, float user_radius) const;
	cube_t get_room_conn_block() const;
	point get_room_conn_pt(float zval) const;
	bool is_blocked_by_gate(point const &p1, point const &p2) const;
	tunnel_seg_t connect_segment_to(point const &pos, bool conn_dir, unsigned new_tseg_ix, bool is_end=0);
	point get_conn_pt(tunnel_conn_t const &c) const;
	cube_t get_conn_bcube(tunnel_conn_t const &c) const;
};
typedef vector<tunnel_seg_t> vect_tunnel_seg_t;

struct breaker_zone_t {
	unsigned rtype, room_start=0, room_end=0;
	int pri_room=-1;
	bool is_dup=0; // used when there are more breakers than rooms, such as in factories with no basements

	breaker_zone_t() : rtype(RTYPE_NOTSET) {} // invalid room
	breaker_zone_t(unsigned t, unsigned s, unsigned e, int pr, bool dup) : rtype(t), room_start(s), room_end(e), pri_room(pr), is_dup(dup) {}
	bool invalid() const {return (rtype != RTYPE_ELEVATOR && room_start == room_end);}
};

struct stairs_landing_base_t : public oriented_cube_t {
	bool bend_dir=0, roof_access=0, stack_conn=0, in_ext_basement=0, in_backrooms=0, against_wall[2]={};
	uint8_t in_mall=0; // 0=not in mall, 1=in mall concourse, 2=in mall back area/hallway
	uint8_t num_stairs=0; // used for malls, since floor spacing may be larger
	stairs_shape shape=SHAPE_STRAIGHT;

	stairs_landing_base_t() {}
	stairs_landing_base_t(cube_t const &c, bool dim_, bool dir_, bool roof_access_, stairs_shape shape_, bool sc=0, bool ieb=0) :
		oriented_cube_t(c, dim_, dir_), roof_access(roof_access_), stack_conn(sc), in_ext_basement(ieb), shape(shape_) {}
	void set_against_wall(bool const val[2]) {against_wall[0] = val[0]; against_wall[1] = val[1];}
	bool is_u_shape        () const {return (shape == SHAPE_U);}
	bool is_l_shape        () const {return (shape == SHAPE_L);}
	bool is_straight       () const {return !(is_u_shape() || is_l_shape());}
	bool has_walled_sides  () const {return (shape == SHAPE_WALLED || shape == SHAPE_WALLED_SIDES);}
	unsigned get_face_id   () const {return (2*dim + dir);}
	unsigned get_num_stairs() const;
	float get_step_length  () const {return get_length()/get_num_stairs();}
	float get_retail_landing_width(float floor_spacing) const {return 0.5*min(get_length(), floor_spacing);}
	float get_stair_dz(float floor_spacing) const {return floor_spacing/(get_num_stairs()+1);}
	float get_wall_hwidth(float floor_spacing) const;
};

struct landing_t : public stairs_landing_base_t {
	bool for_elevator=0, for_ramp=0, has_railing=0, is_at_top=0, not_an_exit=0;
	uint8_t floor_ix=0; // Note: only used for player coll when floor_ix <= 1

	landing_t() {}
	landing_t(cube_t const &c, bool e, uint8_t f, bool dim_, bool dir_,
		bool railing=0, stairs_shape shape_=SHAPE_STRAIGHT, bool roof_access_=0, bool at_top=0, bool sc=0, bool fr=0, bool ieb=0) :
		stairs_landing_base_t(c, dim_, dir_, roof_access_, shape_, sc, ieb), for_elevator(e), for_ramp(fr), has_railing(railing), is_at_top(at_top), floor_ix(f)
	{assert(is_strictly_normalized());}
};

struct stairwell_t : public stairs_landing_base_t {
	bool extends_below=0, extends_above=0;
	uint8_t num_floors=0;
	int16_t stairs_door_ix=-1, not_an_exit_mask=0;

	stairwell_t() {}
	stairwell_t(cube_t const &c, unsigned n, bool dim_, bool dir_, stairs_shape s=SHAPE_STRAIGHT, bool r=0, bool sc=0, bool ieb=0) :
		stairs_landing_base_t(c, dim_, dir_, r, s, sc, ieb), num_floors(n) {}
	stairwell_t(stairs_landing_base_t const &b, unsigned n) : stairs_landing_base_t(b), num_floors(n) {} // can set from landing
};
typedef vector<stairwell_t> vect_stairwell_t;

struct stairs_place_t : public cube_t { // for extended basements
	bool dim, dir, add_railing;
	stairs_place_t(cube_t const &c, bool dim_, bool dir_, bool add_railing_) : cube_t(c), dim(dim_), dir(dir_), add_railing(add_railing_) {}
};

struct escalator_t : public oriented_cube_t { // Note: not yet used
	bool move_dir=0, in_mall=0; // move_dir points upward
	float end_ext=0.0, delta_z=0.0, bot_edge_shift=0.0;

	escalator_t() {}
	escalator_t(cube_t const &c, bool dim_, bool dir_, bool mdir, float ext, float dz, float brz, bool in_mall_=0) :
		oriented_cube_t(c, dim_, dir_), move_dir(mdir), in_mall(in_mall_), end_ext(ext), delta_z(dz), bot_edge_shift(brz) {}
	bool  is_going_up    () const {return move_dir;}
	float get_side_width () const {return 0.10*get_width();}
	float get_side_height() const {return 0.80*get_width();} // ramp/steps to top edge
	float get_upper_hang () const {return 0.25*get_width();}
	float get_floor_thick() const {return 0.01*get_width();}
	float get_side_deltaz() const {return (get_side_height() + bot_edge_shift);} // bottom edge to top edge
	cube_t get_ramp_bcube(bool exclude_sides) const;
	void get_ends_bcube(cube_t &lo_end, cube_t &hi_end, bool exclude_sides) const;
	cube_t get_side_for_end(cube_t const &end, bool lr) const;
	cube_t get_support_pillar() const;
	void get_ramp_bottom_pts(cube_t const &ramp, point bot_pts[4]) const;
	void get_all_cubes(cube_t cubes[7]) const;
};

// door flags
unsigned const DOOR_FLAG_MULT_FLOOR = 0x01; // multi-floor room door
unsigned const DOOR_FLAG_FOR_CLOSET = 0x02; // closet door
unsigned const DOOR_FLAG_BLDG_CONN  = 0x04; // door connecting two different buildings
unsigned const DOOR_FLAG_BACKROOMS  = 0x08; // door in backrooms
unsigned const DOOR_FLAG_SMALL_ROOM = 0x10; // door for a small room such as an elevator equipment room
unsigned const DOOR_FLAG_AUTO_CLOSE = 0x20; // door automatically closes (office bathroom)

struct door_base_t : public cube_t {
	bool dim=0, open_dir=0, hinge_side=0, on_stairs=0;
	uint8_t flags=0;
	uint16_t conn_room[2]={}; // on each side of the door
	// is it useful to store the two rooms in the door/door_stack? this will speed up connectivity searches for navigation and room assignment,
	// but only for finding the second room connected to a door, because we still need to iterate over all doors;
	// unfortunately, it's not easy/cheap to assign these values because the room may not even be added until after the door is placed, so we have to go back and set room1/room2 later
	//uint8_t room1, room2;
	door_base_t() {}
	door_base_t(cube_t const &c, bool dim_, bool dir, bool os=0, bool hs=0) :
		cube_t(c), dim(dim_), open_dir(dir), hinge_side(hs), on_stairs(os) {assert(is_normalized());}
	bool get_check_dirs  () const {return (dim ^ open_dir ^ hinge_side ^ 1);}
	bool get_handle_side () const {return (get_check_dirs() ^ 1);} // 0: handle on left; 1: handle on right
	float get_width      () const {return get_sz_dim(!dim);}
	float get_thickness  () const {return DOOR_THICK_TO_WIDTH*get_width();}
	cube_t get_true_bcube     () const {cube_t bc(*this); bc.expand_in_dim(dim, 0.5*get_thickness()); return bc;}
	cube_t get_clearance_bcube() const {cube_t bc(*this); bc.expand_in_dim(dim,     get_width    ()); return bc;}
	cube_t get_open_door_path_bcube() const;
	cube_t get_open_door_bcube_for_room(cube_t const &room) const;
	bool not_a_room_separator() const {return (on_stairs || get_for_closet() || get_backrooms() || get_bldg_conn() || get_small_room());}
	bool is_same_stack(door_base_t const &d) const {return (d.x1() == x1() && d.y1() == y1());}
	bool is_connected_to_room(unsigned room_id) const {return (!no_room_conn() && (room_id == conn_room[0] || room_id == conn_room[1]) && !get_bldg_conn());}
	bool no_room_conn() const {return (conn_room[0] == 0 && conn_room[1] == 0);}
	bool use_min_open_amt() const {return (get_bldg_conn() || get_small_room());} // building connector doors don't have adj building info to determine open amount
	unsigned get_conn_room(unsigned room_id) const;

	void set_mult_floor() {flags |= DOOR_FLAG_MULT_FLOOR;}
	void set_for_closet() {flags |= DOOR_FLAG_FOR_CLOSET;}
	void set_bldg_conn () {flags |= DOOR_FLAG_BLDG_CONN ;}
	void set_backrooms () {flags |= DOOR_FLAG_BACKROOMS ;}
	void set_small_room() {flags |= DOOR_FLAG_SMALL_ROOM;}
	void set_auto_close() {flags |= DOOR_FLAG_AUTO_CLOSE;}
	void clear_auto_close() {flags &= ~DOOR_FLAG_AUTO_CLOSE;}
	bool get_mult_floor() const {return (flags & DOOR_FLAG_MULT_FLOOR);}
	bool get_for_closet() const {return (flags & DOOR_FLAG_FOR_CLOSET);}
	bool get_bldg_conn () const {return (flags & DOOR_FLAG_BLDG_CONN );}
	bool get_backrooms () const {return (flags & DOOR_FLAG_BACKROOMS );}
	bool get_small_room() const {return (flags & DOOR_FLAG_SMALL_ROOM);}
	bool get_auto_close() const {return (flags & DOOR_FLAG_AUTO_CLOSE);}
};
struct door_stack_t : public door_base_t {
	unsigned first_door_ix=0, num_doors=1; // first_door_ix is on the lowest floor
	door_stack_t() {}
	door_stack_t(door_base_t const &db, unsigned fdix) : door_base_t(db), first_door_ix(fdix) {}
};
struct door_t : public door_base_t {
	bool open=0, blocked=0;
	uint8_t locked=0; // 1=regular lock, >= 2=padlock, where color index is locked-2
	room_type rtype=RTYPE_NOTSET; // room type connected to door, if special; for example a bathroom
	int obj_ix=-1; // for closets, etc.
	float open_amt=0.0; // 0.0=fully closed, 1.0=fully open

	door_t() {}
	door_t(cube_t const &c, bool dim_, bool dir, bool open_=1, bool os=0, bool hs=0) : door_base_t(c, dim_, dir, os, hs), open(open_), open_amt(open ? 1.0 : 0.0) {}
	bool is_closed_and_locked() const {return (!open && locked);}
	bool is_locked_or_blocked(unsigned have_key) const {return (blocked || !check_key_mask_unlocks(have_key));}
	bool is_partially_open() const {return (open_amt != (open ? 1.0 : 0.0));}
	bool is_closet_door   () const {return (obj_ix >= 0 && !is_padlocked());}
	bool is_padlocked     () const {return (locked >= 2);}
	bool is_locked_unlockable() const {return (locked >= MAX_LOCK_INDEX);}
	bool check_key_mask_unlocks(unsigned key_mask) const;
	void set_padlock_color_ix(unsigned ix) {assert(ix < NUM_LOCK_COLORS); locked = ix + 2;}
	void set_locked_unlockable() {locked = MAX_LOCK_INDEX;}; // use a lock for which there is no matching color key
	unsigned get_padlock_color_ix() const {assert(is_padlocked()); assert(locked < MAX_LOCK_INDEX); return (locked - 2);}
	void toggle_open_state(bool allow_partial_open=0);
	void make_auto_close() {set_auto_close(); open = 0; open_amt = 0.0;}
	void make_fully_open_or_closed() {open_amt = (open ? 1.0 : 0.0);}
	bool next_frame();
};
typedef vector<door_stack_t> vect_door_stack_t;
typedef vector<door_t> vect_door_t;

// Note: some of these roof objects are actually on the ground next to houses
enum {ROOF_OBJ_BLOCK=0, ROOF_OBJ_ANT, ROOF_OBJ_WALL, ROOF_OBJ_ECAP, ROOF_OBJ_AC, ROOF_OBJ_SCAP, ROOF_OBJ_SIGN, ROOF_OBJ_SIGN_CONN, ROOF_OBJ_WTOWER, ROOF_OBJ_DUCT,
	ROOF_OBJ_SMOKESTACK, DETAIL_OBJ_COLLIDER, DETAIL_OBJ_COLL_SHAD, DETAIL_OBJ_SHAD_ONLY};
enum {ROOF_TYPE_FLAT=0, ROOF_TYPE_SLOPE, ROOF_TYPE_PEAK, ROOF_TYPE_HIPPED, ROOF_TYPE_DOME, ROOF_TYPE_ONION, ROOF_TYPE_SHED};

struct roof_obj_t : public cube_t {
	uint8_t type;
	roof_obj_t(uint8_t type_=ROOF_OBJ_BLOCK) : type(type_) {}
	roof_obj_t(cube_t const &c, uint8_t type_=ROOF_OBJ_BLOCK) : cube_t(c), type(type_) {assert(is_strictly_normalized());}
};
typedef vector<roof_obj_t> vect_roof_obj_t;

struct indoor_pool_t : cube_t {
	bool valid=0, dim=0, dir=0, bottomless=0;
	int room_ix=-1;
	float shallow_zval=0.0, orig_z1=0.0;
};


// building AI
struct building_loc_t {
	int part_ix=-1, room_ix=-1, stairs_ix=-1; // -1 is not contained; what about elevator_ix?
	unsigned floor_ix=0; // global for this building, rather than the current part/room
	bool operator==(building_loc_t const &loc) const {return (same_room_floor(loc) && part_ix == loc.part_ix && stairs_ix == loc.stairs_ix);}
	bool same_room_floor(building_loc_t const &loc) const {return (room_ix == loc.room_ix && floor_ix == loc.floor_ix);}
};

struct building_dest_t : public building_loc_t {
	int building_ix=-1;
	point pos;
	building_dest_t() {}
	building_dest_t(building_loc_t const &b, point const &pos_, int bix=-1) : building_loc_t(b), building_ix(bix), pos(pos_) {}
	bool is_valid() const {return (building_ix >= 0 && part_ix >= 0 && room_ix >= 0);}
};

struct ext_basement_room_params_t;

struct building_conn_info_t { // use for buildings with connected rooms (for example house extended basements)
	struct conn_room_t : public cube_t {
		unsigned door_ix;
		bool dim, dir, door_is_b;
		conn_room_t(cube_t const &c, unsigned dix, bool dim_, bool dir_, bool dib) : cube_t(c), door_ix(dix), dim(dim_), dir(dir_), door_is_b(dib) {}
	};
	struct conn_pt_t {
		building_t *b;
		vector<conn_room_t> rooms;
		conn_pt_t(building_t *b_) : b(b_) {assert(b != nullptr);}
	};
	vector<conn_pt_t> conn;

	void add_connection(building_t *b, cube_t const &room, unsigned door_ix, bool dim, bool dir, bool door_is_b);
	building_t *get_conn_bldg_for_pt(point const &p, float radius=0.0) const;
	building_t *get_bldg_containing_pt(building_t &parent, point const &p) const;
	bool is_visible_through_conn(building_t const &parent, building_t const &target, vector3d const &xlate, float view_dist, bool expand_for_light=0) const;
	door_t const *get_door_to_conn_part(building_t const &parent, point const &pos_bs) const;
	cube_t get_conn_room_closest_to(building_t const &parent, building_t const &target, point const &pos_bs) const;
	bool point_in_conn_room(point const &pos_bs) const;
};

struct building_walkway_geom_t {
	cube_t bcube;
	bool dim;
	float door_bounds[2][2]={}; // {dir=0, dir=1} x {lo, hi} for exterior doors, one pair per end
	uint8_t has_door=0; // for false doors; one bit per floor

	building_walkway_geom_t(cube_t const &c, bool dim_) : bcube(c), dim(dim_) {}
	bool has_ext_door(bool dir) const {return (door_bounds[dir][0] < door_bounds[dir][1]);}
	float get_length() const {return bcube.get_sz_dim(dim);}
};
struct building_walkway_t : public building_walkway_geom_t { // "owned" walkway, one per connected building
	bool is_owner, open_ends[2]={};
	building_t *conn_bldg;
	vect_cube_with_ix_t windows;
	vect_cube_t frames;
	cube_t bcube_inc_rooms, skyway_conn, elevator_bcube, elevator_cut;

	building_walkway_t(building_walkway_geom_t const &g, bool owner, building_t *b) : building_walkway_geom_t(g), is_owner(owner), conn_bldg(b) {bcube_inc_rooms = bcube;}
	cube_t get_bcube_inc_open_door() const;
	bool has_skyway_conn() const {return !skyway_conn.is_all_zeros();}
	void attach_elevator(cube_t const &e);
};

struct skyway_conn_t : public cube_t {
	bool dim, dir;
	building_t const *const building;
	skyway_conn_t(cube_t const &c, bool dim_, bool dir_, building_t const *b) : cube_t(c), dim(dim_), dir(dir_), building(b) {}
};

struct store_doorway_t : public cube_t {
	unsigned room_id;
	uint8_t move_dir=2; // 0=down, 1=up, 2=stopped
	bool dim, dir, closed, locked=0;
	float open_amt=0.0;

	store_doorway_t(cube_t const &bc, unsigned r, bool dim_, bool dir_, bool c) : cube_t(bc), room_id(r), dim(dim_), dir(dir_), closed(c), open_amt(closed ? 0.0 : 1.0) {}
	bool locked_closed() const {return (closed && locked);}
	float get_gate_z1 () const {return (z1() + open_amt*dz());}
	bool press_button(bool up_dir);
	int next_frame();
};
typedef vector<store_doorway_t> vect_store_doorway_t;

struct ug_elev_info_t {
	cube_t entrance;
	float top_floor_z2;
	bool dim, dir;
	ug_elev_info_t(cube_t const &e, float tfz2, bool dim_, bool dir_) : entrance(e), top_floor_z2(tfz2), dim(dim_), dir(dir_) {}
};
typedef vector<ug_elev_info_t> vect_ug_elev_info_t;

struct pet_tank_t : public cube_t {
	unsigned animal_type, obj_ix;
	pet_tank_t(cube_t const &c, unsigned at, unsigned ix) : cube_t(c), animal_type(at), obj_ix(ix) {}
};

class store_texture_manager_t;

struct building_mall_info_t {
	vect_cube_with_ix_t landings; // ix stores {is_escalator, se_dim, se_dir, ww_dir}
	vect_store_doorway_t store_doorways; // ix stores store room index
	vector<store_info_t> stores; // added during interior object placement
	stairwell_t ent_stairs;
	vect_cube_t bathrooms; // actually bathroom pairs
	vect_cube_t skylights, ext_stairs_elevators;
	vector<pet_tank_t> pet_tanks; // except for fish (rats, snakes, spiders)
	std::unique_ptr<store_texture_manager_t> tmgr;
	cube_t store_bounds, food_court_bounds;
	colorRGBA mall_wall_color=WHITE;
	int city_elevator_ix=-1, ent_stairs_start_ix=-1;
	bool mall_cylin_ducts=0, store_cylin_ducts=0;

	building_mall_info_t();
	~building_mall_info_t();
	void clear_room_details() {stores.clear(); pet_tanks.clear();}
};

struct bldg_industrial_info_t {
	bool entrance_dim, entrance_dir;
	float entrance_pos;
	cube_t floor_space, entrance_area;
	rand_gen_t rgen; // used for generating smoke
	vect_cube_t sub_rooms, pg_extended_pipes;

	// specifically for factories
	struct smoke_source_t {
		point pos;
		float radius=0.0, time=0.0, next_smoke_time=0.0;
		int pid;
		smoke_source_t(point const &pos_, float radius_, int pid_=-1) : pos(pos_), radius(radius_), pid(pid_) {}
	};
	vector2d machine_row_spacing;
	vector<smoke_source_t> smoke_emitters;

	bldg_industrial_info_t(bool dim, bool dir, float epos, cube_t const &fs, cube_t const &ee) :
		entrance_dim(dim), entrance_dir(dir), entrance_pos(epos), floor_space(fs), entrance_area(ee) {}
	void next_frame(particle_manager_t &particle_manager);
	void clear_room_details() {smoke_emitters.clear();}
};

struct building_interior_t {
	vect_cube_t floors, ceilings, fc_occluders, exclusion, open_walls, split_window_walls;
	vect_cube_t walls[2]; // walls are split by dim, which is the separating dimension of the wall
	vect_cube_with_ix_t int_windows; // ix stores room index
	vect_stairwell_t stairwells;
	vect_tunnel_seg_t tunnels;
	vect_door_t doors;
	vect_door_stack_t door_stacks;
	vector<landing_t> landings; // for stairs and elevators
	vector<room_t> rooms;
	vector<elevator_t> elevators;
	vector<escalator_t> escalators;
	vector<person_t> people;
	std::unique_ptr<building_room_geom_t  > room_geom;
	std::unique_ptr<building_nav_graph_t  > nav_graph;
	std::unique_ptr<building_conn_info_t  > conn_info;
	std::unique_ptr<building_mall_info_t  > mall_info;
	std::unique_ptr<bldg_industrial_info_t> ind_info ;
	cube_with_ix_t pg_ramp, attic_access; // ix stores {2*dim + dir}
	indoor_pool_t pool;
	cube_t basement_ext_bcube, elevator_equip_room;
	draw_range_t draw_range;
	unsigned extb_walls_start[2] = {0,0};
	unsigned gen_room_details_pass=0, rgen_seed_ix=0;
	int garage_room=-1, ext_basement_hallway_room_id=-1, ext_basement_door_stack_ix=-1, last_active_door_ix=-1, security_room_ix=-1;
	uint8_t furnace_type=FTYPE_NONE, attic_type=ATTIC_TYPE_RAFTERS;
	bool door_state_updated=0, is_unconnected=0, ignore_ramp_placement=0, placed_people=0, elevators_disabled=0, attic_access_open=0, has_backrooms=0;
	bool elevator_dir=0, extb_wall_dim=0, extb_wall_dir=0, conn_room_in_extb_hallway=0, has_sec_hallways=0;
	uint8_t num_extb_floors=0; // for malls and backrooms
	uint8_t mens_count=0, womens_count=0; // bathrooms
	float water_zval=0.0; // for multilevel backrooms and swimming pools

	building_interior_t();
	~building_interior_t();
	float get_doorway_width() const;
	room_t const &get_room(unsigned room_ix) const {assert(room_ix < rooms.size()); return rooms[room_ix];}
	room_t       &get_room(unsigned room_ix)       {assert(room_ix < rooms.size()); return rooms[room_ix];}
	door_t const &get_door(unsigned door_ix) const {assert(door_ix < doors.size()); return doors[door_ix];}
	door_t       &get_door(unsigned door_ix)       {assert(door_ix < doors.size()); return doors[door_ix];}
	room_t const &get_extb_start_room() const {return get_room(ext_basement_hallway_room_id);} // extb hallway, backrooms, or mall
	room_t       &get_extb_start_room()       {return get_room(ext_basement_hallway_room_id);} // extb hallway, backrooms, or mall
	bool has_mall() const {return bool(mall_info);}
	bool has_mall_ent_stairs() const {return (has_mall() && !mall_info->ent_stairs.is_all_zeros());}
	int get_store_id_for_room(unsigned room_id) const;
	bool is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin=0.0f, bool inc_open=0, bool check_open_dir=0) const;
	bool is_blocked_by_stairs_or_elevator(cube_t const &c, float dmin=0.0f, bool elevators_only=0, int no_check_enter_exit=0) const;
	void get_stairs_and_elevators_bcubes_intersecting_cube(cube_t const &c, vect_cube_t &bcubes, float ends_clearance=0.0, float sides_clearance=0.0) const;
	void sort_for_optimal_culling();
	void remove_excess_capacity();
	void finalize();
	bool update_elevators(building_t const &building, point const &player_pos);
	bool check_sphere_coll(building_t const &building, point &pos, point const &p_last, float radius,
		vect_room_object_t::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix, bool is_ball=0) const;
	bool check_sphere_coll_room_objects(building_t const &building, point &pos, point const &p_last, float radius,
		vect_room_object_t::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix, bool is_ball=0) const;
	room_object_t const &get_elevator_car(elevator_t const &e) const;
	bool check_sphere_coll_walls_elevators_doors(building_t const &building, point &pos, point const &p_last, float radius,
		float wall_test_extra_z, bool is_player, vector3d *cnorm) const;
	bool line_coll(building_t const &building, point const &p1, point const &p2, point &p_int, bool skip_room_geom=0) const;
	point find_closest_pt_on_obj_to_pos(building_t const &building, point const &pos, float pad_dist, bool no_ceil_floor) const;
	void update_dynamic_draw_data() {assert(room_geom); room_geom->update_dynamic_draw_data();}
	void get_avoid_cubes(vect_cube_t &avoid, float z1, float z2, float r_shrink_if_low, float floor_thickness, float floor_ceil_gap,
		bool same_as_player, bool skip_stairs=0, cube_t const *const fires_select_cube=nullptr) const;
	void create_fc_occluders();
	void add_ceil_floor_pair(cube_t cf, float zc, float z, float zf);
	void place_exterior_room(extb_room_t const &room, cube_t const &wall_area, float fc_thick, float wall_thick, ext_basement_room_params_t &P,
		unsigned part_id, unsigned num_lights=0, bool is_hallway=0, unsigned is_building_conn=0, unsigned wall_skip_dim=2, unsigned thin_wall_dir=2);
	colorRGBA get_attic_ceiling_color() const;
	room_t const &get_garage_room() const {assert(garage_room >= 0); return get_room(garage_room);}
	vector<room_t>::const_iterator ext_basement_rooms_start() const;
	bool point_in_U_stairwell(point const &pos, float floor_spacing, bool mall_only=0) const;
	bool point_in_ext_basement_room(point const &pos, float floor_spacing, float expand=0.0) const;
	// tunnels
	bool point_in_tunnel(point const &pos, float expand=0.0) const;
	bool point_near_tunnel_entrance(point const &pos) const;
	int get_tunnel_ix_for_point(point const &pos) const;
	int get_tunnel_ix_for_room (unsigned room_ix) const;
	bool get_tunnel_path_from_room(point const &end_pt,   unsigned  room_ix, float radius, ai_path_t &path) const;
	bool get_tunnel_path_to_room  (point const &start_pt, unsigned &room_ix, ai_path_t &path) const;
	bool get_tunnel_path_two_pts(point const &start_pt, point const &end_pt, float radius, ai_path_t &path) const;
	bool get_tunnel_path(unsigned tix1, unsigned tix2, int prev_tix, ai_path_t &path) const;
	bool cube_in_ext_basement_room(cube_t const &c, bool xy_only) const;
	door_t const &get_ext_basement_door() const;
	void assign_master_bedroom(float window_vspacing, float floor_thickness);
	void assign_door_conn_rooms(unsigned start_ds_ix=0);
	breaker_zone_t get_circuit_breaker_info(unsigned zone_id, unsigned num_zones, float floor_spacing) const;
	void connect_and_add_tunnel_seg(tunnel_seg_t &parent, point const &conn_pt, bool parent_conn_dir, bool child_conn_dir, float gate_dist_from_end);
	void add_tunnel_seg(tunnel_seg_t const &tseg);
}; // end building_interior_t

struct building_stats_t {
	unsigned nbuildings, nparts, ndetails, ntquads, ndoors, ninterior, nrooms, nceils, nfloors, nwalls, nrgeom, nobjs, nverts;
	building_stats_t() : nbuildings(0), nparts(0), ndetails(0), ntquads(0), ndoors(0), ninterior(0), nrooms(0), nceils(0), nfloors(0), nwalls(0), nrgeom(0), nobjs(0), nverts(0) {}
};

struct colored_cube_t;
typedef vector<colored_cube_t> vect_colored_cube_t;
class cube_bvh_t;
class building_indir_light_mgr_t;

struct colored_sphere_t : public sphere_t {
	colorRGBA color;
	colored_sphere_t() {}
	colored_sphere_t(point const &pos_, float radius_, colorRGBA const &color_) : sphere_t(pos_, radius_), color(color_) {}
};

struct ext_step_t : public cube_t {
	bool dim, step_dir, wall_dir, at_door, is_base, at_ground, enclosed, step_up; // enclosed=1 for balconies
	ext_step_t(cube_t const &c, bool dim_, bool sdir, bool wdir, bool door=0, bool base=0, bool ag=0, bool enc=0, bool su=0) :
		cube_t(c), dim(dim_), step_dir(sdir), wall_dir(wdir), at_door(door), is_base(base), at_ground(ag), enclosed(enc), step_up(su) {}
};

struct building_colors_t {
	colorRGBA side_color, wall_color, basement_wall_color, attic_color;
};

struct building_t : public building_geom_t {

	unsigned mat_ix=0;
	uint8_t hallway_dim=2, real_num_parts=0, roof_type=ROOF_TYPE_FLAT; // main hallway dim: 0=x, 1=y, 2=none
	uint8_t roof_dims=0; // for two-part/L-shaped house roofs: 0=auto based on aspect ratio, 1=perpendicular, 2=parallel
	uint8_t street_dir=0; // encoded as 2*dim + dir + 1; 0 is unassigned
	int8_t open_door_ix=-1, basement_part_ix=-1;
	uint8_t has_chimney=0; // 0=none, 1=interior, 2=exterior with fireplace
	uint8_t city_ix=0; // supports up to 256 cities
	uint8_t floor_ext_door_mask=0; // used for multi-family houses
	uint8_t next_unit_id=1; // for apartments and hotels
	uint8_t next_room_num=0; // room numbers per floor
	building_type_t btype=BTYPE_UNSET;
	bool is_house=0, has_garage=0, has_shed=0, has_int_garage=0, has_courtyard=0, has_complex_floorplan=0, has_helipad=0, has_ac=0, has_fake_roof_door=0;
	bool has_tline_conn=0, has_smokestack=0;
	mutable bool has_attic_window=0; // make mutable so that drawing code can update/cache this value
	bool multi_family=0; // apartments, multi-family house, duplex, etc. - split by floor
	bool has_int_fplace=0, has_parking_garage=0, has_small_part=0, has_basement_door=0, has_basement_pipes=0, parts_generated=0, is_in_city=0, has_skylight_light=0;
	bool pri_hall_stairs_to_pg=0, have_walkway_ext_door=0;
	mutable bool has_missing_stairs=0; // only used for printing a warning
	uint8_t retail_floor_levels=0;
	int8_t courtyard_door_ix=-1;
	mutable bool player_visited=0; // for stats tracking
	colorRGBA side_color=WHITE, roof_color=WHITE, detail_color=BLACK, door_color=WHITE, wall_color=WHITE;
	cube_t bcube, coll_bcube, pri_hall, driveway, porch, assigned_plot, exterior_flag, ladder, deck_bounds;
	mutable cube_t city_driveway; // set by city gen, which only has a const ref to the building; technically this is cached city state, and not directly used by the building
	vect_cube_t parts, fences;
	vect_cube_with_ix_t skylights, gutters;
	vect_roof_obj_t details; // cubes on the roof - antennas, AC units, etc.
	vect_tquad_with_ix_t roof_tquads, doors;
	vector<colored_sphere_t> ext_lights;
	vector<vect_point> per_part_ext_verts; // only used for non-cube buildings
	vector<ext_step_t> ext_steps;
	vector<building_walkway_t> walkways;
	std::shared_ptr<building_interior_t> interior;
	std::string name; // company name for office building; family name for house
	std::string address; // only used for city buildings on roads
	vertex_range_t ext_side_qv_range;
	point tree_pos; // (0,0,0) is unplaced/no tree
	float ao_bcz2=0.0, ground_floor_z1=0.0, interior_z2=0.0, water_damage=0.0, crack_damage=0.0;

	friend class building_indir_light_mgr_t;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_) {}
	building_t(building_geom_t const &bg) : building_geom_t(bg) {}
	static float get_scaled_player_radius();
	static float get_min_front_clearance () {return 2.05f*get_scaled_player_radius();} // slightly larger than the player diameter
	static float get_min_front_clearance_inc_people();
	bool is_valid() const {return !bcube.is_all_zeros();}
	bool has_interior () const {return bool(interior);}
	bool has_conn_info() const {return (interior && interior->conn_info);}
	bool has_room_geom() const {return (has_interior() && interior->room_geom);}
	bool has_sec_bldg () const {return (has_garage || has_shed);}
	bool has_pri_hall () const {return (hallway_dim <= 1);} // otherwise == 2
	bool has_basement () const {return (basement_part_ix >= 0);}
	bool has_driveway () const {return !driveway.is_all_zeros();}
	bool has_a_garage () const {return (has_garage || has_int_garage);} // external or internal
	bool has_attic    () const {return (interior && !interior->attic_access.is_all_zeros());}
	bool has_porch    () const {return !porch.is_all_zeros();}
	bool has_people   () const {return (interior && !interior->people.empty());}
	bool has_retail   () const {return (retail_floor_levels > 0);}
	bool has_tall_retail() const {return (retail_floor_levels > 1);}
	bool has_mall       () const {return (interior && interior->has_mall());}
	bool has_mall_ent_stairs() const {return (interior && interior->has_mall_ent_stairs());}
	bool has_mall_skylight  () const {return (has_mall() && !interior->mall_info->skylights.empty());}
	bool is_office_bldg () const {return (btype == BTYPE_OFFICE    );}
	bool is_apartment   () const {return (btype == BTYPE_APARTMENT );}
	bool is_hotel       () const {return (btype == BTYPE_HOTEL     );}
	bool is_factory     () const {return (btype == BTYPE_FACTORY   );}
	bool is_hospital    () const {return (btype == BTYPE_HOSPITAL  );}
	bool is_school      () const {return (btype == BTYPE_SCHOOL    );}
	bool is_parking     () const {return (btype == BTYPE_PARKING   );}
	bool is_warehouse   () const {return (btype == BTYPE_WAREHOUSE );}
	bool is_powerplant  () const {return (btype == BTYPE_POWERPLANT);}
	bool is_police_stat () const {return (btype == BTYPE_POLICE    );}
	bool is_fire_stat   () const {return (btype == BTYPE_FIRE_STAT );}
	bool is_apt_or_hotel() const {return (is_apartment() || is_hotel());}
	bool is_residential () const {return (is_house || is_apt_or_hotel());}
	bool is_industrial  () const {return (is_factory() || is_warehouse() || is_powerplant());}
	bool is_commercial  () const {return (!is_residential() && !is_industrial());}
	bool is_retail_part(cube_t const &part) const {return (has_retail() && part.z1() == ground_floor_z1);}
	bool skip_top_of_ceilings() const {return (roof_type == ROOF_TYPE_FLAT || !is_house || has_attic());}
	bool enable_driveway_coll() const {return !is_rotated();} // no collision with rotated driveways/porches for now
	bool has_pg_ramp() const {return (interior && !interior->pg_ramp.is_all_zeros());}
	bool can_extend_stairs_to_pg(unsigned &stairs_ix) const;
	bool is_basement(vect_cube_t::const_iterator it) const {return (int(it - parts.begin()) == basement_part_ix);}
	bool is_pos_in_basement(point const &pos) const {return ((has_basement() && parts[basement_part_ix].contains_pt(pos)) || point_in_extended_basement(pos));}
	bool maybe_has_ext_door_this_floor(float part_z1, unsigned floor_ix) const;
	void get_garage_dim_dir(cube_t const &garage, bool &dim, bool &dir) const;
	unsigned get_attic_part_ix   () const;
	room_t const &get_retail_room() const {assert(interior && !interior->rooms.empty()); assert(has_retail()); return interior->rooms.front();} // always the first room
	cube_t const &get_retail_part() const {assert(has_retail()); assert(!parts.empty()); return parts.front();} // always the first part
	cube_t const &get_basement   () const {assert(has_basement()); return parts[basement_part_ix   ];}
	cube_t const &get_attic_part () const {assert(has_attic   ()); return parts[get_attic_part_ix()];}
	bool get_retail_long_dim     () const;
	int check_player_in_basement(point const &pos) const;
	colorRGBA get_avg_side_color  () const {return side_color  .modulate_with(get_material().side_tex.get_avg_color());}
	colorRGBA get_avg_roof_color  () const {return roof_color  .modulate_with(get_material().roof_tex.get_avg_color());}
	colorRGBA get_avg_detail_color() const {return detail_color.modulate_with(get_material().roof_tex.get_avg_color());}
	colorRGBA get_door_handle_color() const;
	building_mat_t const &get_material() const;
	bool has_windows         () const {return get_material().add_windows;}
	bool has_int_windows     () const {return (DRAW_CITY_INT_WINDOWS || has_windows());}
	float get_floor_thick_val() const {return (is_house ? FLOOR_THICK_VAL_HOUSE : (has_windows() ? FLOOR_THICK_VAL_OFFICE : FLOOR_THICK_VAL_WINDOWLESS));}
	float get_elevator_fc_thick_scale() const {return 1.005*0.5*get_floor_thick_val();}
	float get_window_vspace  () const {return get_material().get_floor_spacing();}
	float get_floor_thickness() const {return get_floor_thick_val()*get_window_vspace();}
	float get_fc_thickness   () const {return 0.5*get_floor_thickness();} // floor/ceiling thickness
	float get_wall_thickness () const {return WALL_THICK_VAL*get_window_vspace();}
	float get_wind_trim_thick() const {return 0.75*get_wall_thickness();}
	float get_trim_thickness () const {return 0.1 *get_wall_thickness();}
	float get_trim_height    () const {return 0.04*get_window_vspace ();}
	float get_floor_ceil_gap () const {return (get_window_vspace() - get_floor_thickness());}
	float get_door_height    () const {return 0.95f*get_floor_ceil_gap();} // set height based on window spacing, 95% of ceiling height (may be too large)
	float get_office_bldg_door_height() const {return 1.06*get_door_height();} // a bit taller
	float get_attic_beam_depth()const {return 0.08*get_window_vspace();}
	float get_min_wall_len   () const {return 2.0 *get_window_vspace();}
	float get_door_shift_dist() const {return 0.01*get_window_vspace();}
	float get_rug_thickness  () const {return 0.0010*get_window_vspace();}
	float get_flooring_thick () const {return 0.0012*get_window_vspace();}
	float get_doorway_width  () const;
	float get_landing_width  () const {return 1.0*get_doorway_width();} // for L-shaped stairs
	float get_nominal_doorway_width   () const {return DOOR_WIDTH_SCALE*get_window_vspace();} // constant per-building, but not exactly the same as get_doorway_width()
	float get_office_ext_doorway_width() const {return DOOR_WIDTH_SCALE_OFFICE*get_window_vspace();}
	float get_industrial_window_z1    () const;
	bool is_ground_floor_excluding_retail(float zval) const;
	float get_ground_floor_z_thresh(bool for_spider) const;
	void gen_rotation(rand_gen_t &rgen);
	point get_inv_rot_pos(point const &pos) const;
	void maybe_inv_rotate_pos_dir(point &pos, vector3d &dir) const;
	void set_z_range(float z1, float z2);
	bool check_part_contains_pt_xy(cube_t const &part, unsigned part_id, point const &pt) const;
	bool check_part_contains_cube_xy(cube_t const &part, unsigned part_id, cube_t const &c) const;
	bool cube_int_parts_no_sec(cube_t const &c) const;
	bool check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs) const;
	bool check_cube_within_part_sides(cube_t const &c) const;
	bool check_pt_within_part_sides(point const &p) const;
	bool check_pt_in_retail_room(point const &p) const;
	bool check_pt_in_or_near_walkway(point const &p, bool owned_only, bool inc_open_door, bool inc_conn_room) const;
	bool is_connected_with_walkway(building_t const &target, float zval=0.0) const;
	vect_cube_t::const_iterator get_real_parts_end() const {return (parts.begin() + real_num_parts);}
	vect_cube_t::const_iterator get_real_parts_end_inc_sec() const {return (get_real_parts_end() + has_sec_bldg());}
	vect_point const &get_part_ext_verts(unsigned part_id) const {assert(part_id < per_part_ext_verts.size()); return per_part_ext_verts[part_id];}
	cube_t const &get_sec_bldg () const {assert(has_sec_bldg()); assert(real_num_parts < parts.size()); return parts[real_num_parts];}
	cube_t const &get_chimney  () const {assert(has_chimney      && parts.size() > 1); return parts.back();}
	cube_t const &get_fireplace() const {assert(has_chimney == 2 && parts.size() > 2); return parts[parts.size()-2];}
	cube_t get_interior_bcube(bool inc_ext_basement) const;
	cube_t get_ext_vis_bcube() const;
	void union_with_coll_bcube(cube_t const &c);

	bool check_sphere_coll(point const &pos, float radius, bool xy_only, vector3d *cnorm=nullptr) const {
		point pos2(pos);
		return check_sphere_coll(pos2, pos, zero_vector, radius, xy_only, cnorm);
	}
	bool check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector3d *cnorm=nullptr, bool check_interior=0) const;
	bool check_sphere_coll_inner(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector3d *cnorm=nullptr, bool check_interior=0) const;
	bool check_sphere_coll_interior(point &pos, point const &p_last, float radius, bool is_in_attic, bool xy_only, vector3d *cnorm) const;
	bool check_cube_intersect_non_main_part(cube_t const &c) const;
	bool point_in_elevator   (point const &pos, bool check_elevator_car=0) const;
	bool point_in_escalator  (point const &pos) const {return (interior && check_vect_cube_contains_pt(interior->escalators, pos));}
	bool point_in_stairwell  (point const &pos) const {return (interior && check_vect_cube_contains_pt(interior->stairwells, pos));}
	bool point_in_U_stairwell(point const &pos) const {return (interior && interior->point_in_U_stairwell(pos, get_window_vspace()));}
	bool check_pos_in_unlit_room(point const &pos) const;
	bool check_pos_in_unlit_room_recur(point const &pos, std::set<unsigned> &rooms_visited, int known_room_id=-1) const;
	bool is_room_windowless(room_t const &room) const;
	bool are_rooms_connected(room_t const &r1, room_t const &r2, float zval, bool check_door_open) const;
	bool all_room_int_doors_closed(unsigned room_ix, float zval) const;
	unsigned check_line_coll(point const &p1, point const &p2, float &t, bool occlusion_only=0, bool ret_any_pt=0, bool no_coll_pt=0, bool check_non_coll=0) const;
	bool get_interior_color_at_xy(point const &pos, colorRGBA &color) const;
	bool point_in_mall_elevator_entrance(point const &pos, bool inc_front_space) const;
	int check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points,
		bool inc_attic=0, bool inc_ext_basement=0, bool inc_roof_acc=0, bool inc_details=0, cube_t *coll_cube=nullptr) const;
	bool point_under_attic_roof(point const &pos, vector3d *const cnorm=nullptr) const;
	bool point_in_attic        (point const &pos, vector3d *const cnorm=nullptr) const;
	bool cube_in_attic         (cube_t const &c) const;
	bool point_in_mall         (point const &pos) const {return (has_mall     () && point_in_extended_basement_not_basement(pos));} // including stores and hallways
	bool point_in_industrial   (point const &pos) const {return (is_industrial() && get_industrial_area().contains_pt(pos));}
	bool check_point_xy_in_part(point const &pos) const;
	bool player_can_see_outside() const;
	bool player_can_see_mall_skylight(vector3d const &xlate) const;
	bool player_can_see_out_mall_skylight(vector3d const &xlate) const;
	bool player_can_see_in_mall_skylight(vector3d const &xlate) const;
	void set_building_colors(building_colors_t &bcolors) const;
	bool ray_cast_exterior_walls(point const &p1, point const &p2, vector3d &cnorm, float &t) const;
	bool ray_cast_interior(point const &pos, vector3d const &dir, cube_t const &valid_area, cube_bvh_t const &bvh, bool in_attic, bool in_ext_basement,
		building_colors_t const &bcolors, point &cpos, vector3d &cnorm, colorRGBA &ccolor, rand_gen_t *rgen=nullptr) const;
	void create_building_volume_light_texture(unsigned bix, point const &target, unsigned &tid) const;
	bool ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const;
	cube_t calc_parts_bcube() const;
	cube_t get_unrotated_parts_bcube() const {return (is_rotated() ? calc_parts_bcube() : bcube);}
	void calc_bcube_from_parts();
	void adjust_part_zvals_for_floor_spacing(cube_t &c) const;
	void gen_geometry(int rseed1, int rseed2);
	void setup_damage_vals();
	cube_t place_door(cube_t const &base, bool dim, bool dir, float door_height, float door_center, float door_pos,
		float door_center_shift, float width_scale, bool can_fail, bool opens_up, rand_gen_t &rgen, unsigned floor_ix=0) const;
	bool check_walkway_door_clearance(cube_t const &c, bool dim) const;
	bool add_walkway_door(building_walkway_geom_t &walkway, bool dir, unsigned part_ix);
	void gen_house(cube_t const &base, rand_gen_t &rgen);
	bool maybe_add_house_driveway(cube_t const &plot, unsigned building_ix) const;
	bool get_power_point(vector<point> &ppts) const;
	void add_solar_panels(rand_gen_t &rgen);
	bool add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_office_building, bool roof_access=0, bool courtyard=0, bool for_walkway=0);
	float gen_peaked_roof(cube_t const &top_, float peak_height, bool dim, float extend_to, float max_dz, unsigned skip_side_tri);
	float gen_hipped_roof(cube_t const &top_, float peak_height, float extend_to);
	float gen_sloped_roof_for_stacked_parts(cube_t const &bot, cube_t const &top);
	void place_roof_ac_units(unsigned num, float sz_scale, cube_t const &bounds, vect_cube_t const &avoid, rand_gen_t &rgen);
	void add_roof_walls(cube_t const &c, float wall_width, bool overlap_corners, cube_t out[4]);
	void gen_details(rand_gen_t &rgen, bool is_rectangle);
	void add_smokestack(rand_gen_t &rgen);
	void maybe_add_skylight(rand_gen_t &rgen);
	void add_company_sign(rand_gen_t &rgen);
	cube_t get_helipad_bcube() const;
	int get_num_windows_on_side(float xy1, float xy2) const;
	int get_num_windows_on_side(cube_t const &c, bool dim) const {return get_num_windows_on_side(c.d[dim][0], c.d[dim][1]);}
	float get_window_h_border() const;
	float get_window_v_border() const;
	float get_hspacing_for_part(cube_t const &part, bool dim) const;
	bool interior_enabled() const;
	void gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes);
	void assign_special_room_types(vector<unsigned> &utility_room_cands, vector<unsigned> &special_room_cands, unsigned doors_start, rand_gen_t &rgen);
	void add_conference_room_window(unsigned room_ix);
	void divide_last_room_into_apt_or_hotel(unsigned room_row_ix, unsigned hall_num_rooms, unsigned tot_num_windows,
		unsigned windows_per_room, unsigned windows_per_room_side, bool hall_dim, bool hall_dir, rand_gen_t &rgen);
	void add_hospital_bathrooms(unsigned rooms_start, rand_gen_t &rgen);
	bool maybe_create_nested_bathroom(room_t &room, rand_gen_t &rgen);
	void create_industrial_floorplan(unsigned part_id, float window_hspacing[2], float window_border, rand_gen_t &rgen);
	bool maybe_assign_interior_garage(bool &gdim, bool &gdir);
	void add_parking_garage_ramp(rand_gen_t &rgen);
	bool add_machines_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool less_clearance=0);
	void add_machines_to_factory(rand_gen_t rgen, room_t const &room, cube_t const &place_area, float zval,
		unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned objs_start_inc_beams, cube_t const &ladder);
	void add_chem_tank_gauge(room_object_t const &tank, float radius, float height);
	void add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned part_ix, unsigned num_floors,
		unsigned rooms_start, bool use_hallway, bool first_part_this_stack, float window_hspacing[2], float window_border, bool is_single_floor);
	void connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part, unsigned lower_part_ix);
	void create_two_story_tall_rooms(rand_gen_t &rgen);
	void setup_courtyard();
	void gen_room_details(rand_gen_t &rgen, unsigned building_ix);
	unsigned calc_floor_offset(float zval, float floor_spacing) const;
	cube_t get_init_elevator_car(elevator_t const &elevator) const;
	void add_stairs_and_elevators(rand_gen_t &rgen);
	int get_ext_door_dir(cube_t const &door_bcube, bool dim) const;
	bool add_sign_by_door(tquad_with_ix_t const &door, bool outside, std::string const &text, colorRGBA const &color, bool emissive);
	void add_doorbell_lamp_and_porch_items(tquad_with_ix_t const &door, rand_gen_t &rgen);
	void add_exterior_door_items(rand_gen_t &rgen);
	unsigned get_street_house_number() const;
	void gen_building_doors_if_needed(rand_gen_t &rgen);
	void maybe_add_special_roof(rand_gen_t &rgen);
	void gen_sloped_roof(rand_gen_t &rgen, cube_t const &top);
	void add_roof_to_bcube();
	void gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax);
	tid_nm_pair_t get_basement_wall_texture  () const;
	tid_nm_pair_t get_industrial_wall_texture() const;
	tid_nm_pair_t get_attic_texture() const;
	tid_nm_pair_t get_interior_ext_wall_texture() const {return (is_industrial() ? get_industrial_wall_texture() : get_material().wall_tex);}
	tid_nm_pair_t get_tile_floor_texture       () const;
	colorRGBA get_floor_tex_and_color(cube_t const &floor_cube, tid_nm_pair_t &tex) const;
	colorRGBA get_ceil_tex_and_color (cube_t const &ceil_cube,  tid_nm_pair_t &tex) const;
	colorRGBA get_trim_color() const {return (is_house ? WHITE : DK_GRAY);}
	bool has_tile_floor() const;
	void get_all_drawn_exterior_verts(building_draw_t &bdraw);
	void get_detail_shadow_casters   (building_draw_t &bdraw);
	void get_all_drawn_ext_wall_verts(building_draw_t &bdraw);
	void get_basement_ext_wall_verts (building_draw_t &bdraw) const;
	void get_all_drawn_interior_verts(building_draw_t &bdraw);
	void get_walkway_interior_verts  (building_draw_t &bdraw, building_walkway_t &w);
	void get_all_drawn_window_verts  (building_draw_t &bdraw, bool lights_pass=0, float offset_scale=1.0,
		point const *only_cont_pt_in=nullptr, bool no_skylights=0, bool draw_int_windows=0) const;
	void get_all_drawn_window_verts_as_quads(vect_vnctcc_t &verts) const;
	bool get_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, point const &pos, vector3d const &view_dir, float dist, bool update_state, bool only_open);
	void get_ext_door_verts(building_draw_t &bdraw, point const &viewer, vector3d const &view_dir, int skip_door_ix) const;
	bool get_all_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, vector<point> const &pts, float dist);
	void player_not_near_building() {register_open_ext_door_state(-1);}
	int find_ext_door_close_to_point(tquad_with_ix_t &door, point const &pos, float dist) const;
	bool point_near_ext_door(point const &pos, float dist) const;
	bool get_building_door_pos_closest_to(point const &target_pos, point &door_pos, bool inc_garage_door) const;
	cube_t register_deck_and_get_part_bounds(cube_t const &deck);
	void get_split_int_window_wall_verts(building_draw_t &bdraw_front, building_draw_t &bdraw_back, point const &only_cont_pt_in, bool make_all_front=0) const;
	void get_ext_wall_verts_no_sec(building_draw_t &bdraw) const;
	void get_walkway_end_verts(building_draw_t &bdraw, point const &pos) const;
	void write_basement_entrance_depth_pass(shader_t &s) const;
	void add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, bool sec_camera_mode,
		occlusion_checker_noncity_t const &oc, vect_cube_with_ix_t &ped_bcubes, cube_t &lights_bcube);
	colorRGBA get_retail_light_color() const;
	void run_light_motion_detect_logic(point const &camera_bs);
	bool toggle_room_light(point const &closest_to, bool sound_from_closest_to=0, int room_id=-1, bool inc_lamps=1, bool closet_light=0, bool known_in_attic=0);
	bool toggle_walkway_lights(point const &pos);
	void toggle_light_object(room_object_t const &light, point const &sound_pos);
	void register_light_state_change(room_object_t const &light, point const &sound_pos, bool is_lamp=0);
	void toggle_circuit_breaker(bool is_on, unsigned zone_id, unsigned num_zones);
	bool chair_can_be_rotated(room_object_t const &chair) const;
	void run_player_interact_logic(point const &camera_bs);
	bool apply_player_action_key(point const &closest_to_in, vector3d const &in_dir_in, int mode, bool check_only=0, bool no_check_conn_building=0);
	void assign_correct_room_to_object(room_object_t &obj) const;
	bool move_nearest_object(point const &at_pos, vector3d const &in_dir, float range, int mode);
	bool interact_with_object(unsigned obj_ix, point const &int_pos, point const &query_ray_end, vector3d const &int_dir);
	bool adjust_blinds_state(unsigned obj_ix);
	void add_box_contents(room_object_t const &box);
	void toggle_door_state(unsigned door_ix, bool player_in_this_building, bool by_player, point const &actor_pos);
	void notify_door_fully_closed_state(door_t const &door);
	void handle_items_intersecting_closed_door(door_t const &door);
	void doors_next_frame(point const &player_pos);
	bool set_room_light_state_to(room_t const &room, float zval, bool make_on);
	void set_obj_lit_state_to(unsigned room_id, float light_z2, bool lit_state);
	bool player_pickup_object(point const &at_pos, vector3d const &in_dir);
	void register_player_change_floor(unsigned old_floor, unsigned new_floor) const;
	void register_player_enter_building() const;
	void register_player_exit_building (bool entered_another_building) const;
	bool check_for_wall_ceil_floor_int(point const &p1, point const &p2, bool inc_pg_br_walls=1) const;
	bool line_intersect_stairs_or_ramp(point const &p1, point const &p2) const;
	bool check_cube_on_or_near_stairs(cube_t const &c) const;
	bool drop_room_object(room_object_t &obj, point const &dest, point const &player_pos, bool dim, bool dir);
	bool maybe_use_last_pickup_room_object(point const &player_pos, bool no_time_check=0, bool random_dir=0);
	bool maybe_update_tape(point const &player_pos, bool end_of_tape);
	void handle_vert_cylin_tape_collision(point &cur_pos, point const &prev_pos, float z1, float z2, float radius, bool is_player) const;
	void draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, occlusion_checker_noncity_t &oc, vector3d const &xlate,
		unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building);
	void gen_and_draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, shader_t &amask_shader, occlusion_checker_noncity_t &oc, vector3d const &xlate, unsigned building_ix,
		bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building, bool ext_basement_conn_visible, bool mall_visible);
	bool has_glass_floor() const {return (has_room_geom() && !interior->room_geom->glass_floors.empty());}
	bool glass_floor_visible(vector3d const &xlate, bool from_outside_building=0) const;
	bool point_over_glass_floor(point const &pos, bool inc_escalator=0) const;
	void draw_glass_surfaces(vector3d const &xlate) const;
	void draw_factory_alpha (vector3d const &xlate) const;
	bool has_cars_to_draw(bool player_in_building) const;
	void draw_cars_in_building(shader_t &s, vector3d const &xlate, bool player_in_this_building, bool shadow_only) const;
	bool check_for_water_splash(point const &pos_bs, float size=1.0, bool full_room_height=0, bool draw_splash=0, bool alert_zombies=1) const;
	cube_t calc_splash_bounds(point const &pos) const;
	void draw_water(vector3d const &xlate) const;
	void debug_people_in_building(shader_t &s, point const &camera_bs) const;
	void subtract_stairs_and_elevators_from_cube(cube_t const &c, vect_cube_t &cube_parts, bool inc_stairs=1, bool inc_elevators=1) const;
	void add_split_roof_shadow_quads(building_draw_t &bdraw) const;
	void clear_room_geom(bool even_if_player_modified=0);
	void clear_and_regen_new_seed();
	void update_grass_exclude_at_pos(point const &pos, vector3d const &xlate, bool camera_in_building) const;
	void add_signs(vector<sign_t> &signs) const;
	void add_flags(vector<city_flag_t> &flags);
	void update_stats(building_stats_t &s) const;
	bool are_rooms_connected_without_using_room_or_door(unsigned room1, unsigned room2, unsigned room_exclude, int door_exclude=-1, float door_ex_zval=0.0) const;
	bool is_room_adjacent_to_ext_door(cube_t const &room, bool front_door_only=0) const;
	bool cube_int_ext_door(cube_t const &c) const;
	room_t const &get_room(unsigned room_ix) const {assert(interior); return interior->get_room(room_ix);}
	room_t       &get_room(unsigned room_ix)       {assert(interior); return interior->get_room(room_ix);}
	door_t const &get_door(unsigned door_ix) const {assert(interior); return interior->get_door(door_ix);}
	door_t       &get_door(unsigned door_ix)       {assert(interior); return interior->get_door(door_ix);}
	point get_center_of_room(unsigned room_ix) const {return get_room(room_ix).get_cube_center();}
	room_t const &get_pool_room() const {assert(interior); return get_room(interior->pool.room_ix);}
	std::string get_room_name(point const &pos, int room_id=-1, unsigned floor_ix=0) const;
	void create_pending_textures() const;
	void clear_clothing_textures() const;
	bool bind_custom_clothing_texure(room_object_t const &obj) const;

	// building AI people
	unsigned count_connected_room_components();
	bool place_people_if_needed(unsigned building_ix, float radius, vect_point &locs) const;
	void all_ai_room_update(rand_gen_t &rgen, float delta_dir);
	int ai_room_update(person_t &person, float delta_dir, unsigned person_ix, rand_gen_t &rgen);
	int run_ai_elevator_logic(person_t &person, float delta_dir, rand_gen_t &rgen);
	bool run_ai_pool_logic  (person_t &person, float &speed_mult) const;
	bool run_ai_tunnel_logic(person_t &person, float &speed_mult) const;
	bool maybe_zombie_retreat(unsigned person_ix, point const &hit_pos);
	void register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity);
	bool is_room_backrooms(unsigned room_ix) const {return get_room(room_ix).is_backrooms();}
	bool is_single_large_room         (int room_ix) const {return (room_ix >= 0 && get_room(room_ix).is_single_large_room());}
	bool is_single_large_room_or_store(int room_ix) const {return (room_ix >= 0 && get_room(room_ix).is_single_large_room_or_store());}
	bool is_above_retail_area(point const &pos) const;
	bool is_pos_in_pg_or_backrooms(point const &pos) const;
	bool has_backrooms_or_mall() const {return (interior && (interior->has_backrooms || has_mall()));}
	point get_retail_upper_stairs_landing_center() const;
	void add_sub_rooms_to_avoid_if_needed(room_t const &room, vect_cube_t &avoid) const;
private:
	void clear_existing_room_geom();
	void build_nav_graph() const;
	bool is_valid_ai_placement(point const &pos, float radius, bool skip_nocoll, bool no_check_objs=0) const;
	bool choose_dest_goal(person_t &person, rand_gen_t &rgen) const;
	int  choose_dest_room(person_t &person, rand_gen_t &rgen) const;
	int  maybe_use_escalator(person_t &person, building_loc_t const &loc, bool last_used_escalator, rand_gen_t &rgen) const;
	bool select_person_dest_in_room(person_t &person, rand_gen_t &rgen, room_t const &room) const;
	void get_avoid_cubes(float zval, float height, float radius, vect_cube_t &avoid, bool following_player, cube_t const *const fires_select_cube=nullptr) const;
	bool find_route_to_point(person_t &person, float radius, bool is_first_path, bool following_player, ai_path_t &path) const;
	void add_escalator_points(person_t const &person, ai_path_t &path) const;
	bool stairs_contained_in_part(stairwell_t const &s, cube_t const &p) const;
	bool no_stairs_exit_on_floor(stairwell_t const &stairs, float zval) const;
	void find_nearest_stairs_ramp_esc(point const &p1, point const &p2, vector<unsigned> &nearest_stairs, int part_ix=-1) const;
	int find_nearest_elevator_this_floor(point const &pos) const;
	void ai_room_lights_update(person_t const &person);
	void move_person_to_not_collide(person_t &person, person_t const &other, point const &new_pos, float rsum, float coll_dist) const;
	void register_player_hiding(room_object_t const &hiding_obj) const;
	elevator_t       &get_elevator(unsigned eix)       {assert(interior); assert(eix < interior->elevators.size()); return interior->elevators[eix];}
	elevator_t const &get_elevator(unsigned eix) const {assert(interior); assert(eix < interior->elevators.size()); return interior->elevators[eix];}

	// animals
public:
	template<typename T> void add_animals_on_floor(T &animals, unsigned building_ix, unsigned num_min, unsigned num_max, float sz_min, float sz_max) const;
	void update_animals      (point const &camera_bs, unsigned building_ix);
	void update_rats         (point const &camera_bs, unsigned building_ix);
	void update_sewer_rats   (point const &camera_bs, unsigned building_ix);
	void update_pet_rats     (point const &camera_bs, unsigned building_ix);
	void update_spiders      (point const &camera_bs, unsigned building_ix);
	void update_sewer_spiders(point const &camera_bs, unsigned building_ix);
	void update_snakes       (point const &camera_bs, unsigned building_ix);
	void update_pet_snakes   (point const &camera_bs, unsigned building_ix);
	void update_pet_birds    (point const &camera_bs, unsigned building_ix);
	void update_insects      (point const &camera_bs, unsigned building_ix);
	void get_objs_at_or_below_ground_floor(vect_room_object_t &ret, bool for_spider) const;
	bool begin_fish_draw() const;
	void rat_bite_player(point const &pos, float damage, rand_gen_t &rgen);
private:
	// animals
	point gen_animal_floor_pos(float radius, bool place_in_attic, bool not_player_visible, bool pref_dark_room, bool not_by_ext_door, rand_gen_t &rgen) const;
	bool add_rat(point const &pos, float hlength, vector3d const &dir, point const &placed_from, bool &dead);
	void update_rat(rat_t &rat, point const &camera_bs, float timestep, float &max_xmove, bool can_attack_player, rand_gen_t &rgen);
	void scare_rat(rat_t &rat, point const &camera_bs) const;
	void scare_rat_at_pos(rat_t &rat, point const &scare_pos, float amount, bool by_sight) const;
	bool update_spider_pos_orient(spider_t &spider, point const &camera_bs, float timestep, rand_gen_t &rgen) const;
	void update_spider(spider_t &spider, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen);
	bool maybe_squish_animals(room_object_t const &obj, point const &player_pos);
	void update_snake(snake_t  &snake,  point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen);
	int  check_for_animal_coll(building_animal_t const &A, float hheight, float z_center_offset, bool on_floor_only, bool skip_player,
		point const &camera_bs, float timestep, point const &old_pos, point const &query_pos, vector3d &coll_dir) const;
	int  check_for_snake_coll(snake_t const &snake, point const &camera_bs, float timestep, point const &old_pos, point const &query_pos, vector3d &coll_dir) const;
	void update_insect(insect_t &insect, point const &camera_bs, float timestep, vector<pair<float, point>> &targets, rand_gen_t &rgen) const;
	void update_fly   (insect_t &fly,    point const &camera_bs, float timestep, vector<pair<float, point>> &targets, rand_gen_t &rgen) const;
	void update_roach (insect_t &roach,  point const &camera_bs, float timestep, rand_gen_t &rgen) const;
	void maybe_bite_and_poison_player(point const &pos, point const &camera_bs, vector3d const &dir, float coll_radius, float damage, int poison_type, rand_gen_t &rgen);

	bool is_pos_inside_building(point const &pos, float xy_pad, float hheight, bool inc_attic=1) const;
	void get_room_obj_cubes(room_object_t const &c, point const &pos, vect_cube_t &lg_cubes, vect_cube_t &sm_cubes, vect_cube_t &non_cubes) const;
	int  check_line_coll_expand(point const &p1, point const &p2, float radius, float hheight, bool for_spider=0) const;
	bool check_line_of_sight_large_objs(point const &p1, point const &p2) const;
	bool check_line_int_interior_window(point const &p1, point const &p2) const;
	bool check_and_handle_dynamic_obj_coll(point &pos, point const &cur_obj_pos, float radius,
		float z1, float z2, point const &camera_bs, bool for_spider, bool skip_player=0) const;
	bool get_begin_end_room_objs_on_ground_floor(float zval, bool for_spider, vect_room_object_t::const_iterator &b, vect_room_object_t::const_iterator &e) const;
public:
	int get_room_containing_pt(point const &pt) const;
	int get_room_containing_camera(point const &camera_rot) const;
	unsigned get_attic_room_id() const;
	int room_or_adj_room_has_stairs(int room_ix, float zval, bool inc_adj_rooms, bool check_door_open) const;
	void register_player_in_building(point const &camera_bs, unsigned building_id) const;
	bool maybe_teleport_to_screenshot() const;
	void update_security_cameras(point const &camera_bs);
	bool check_if_against_window(cube_t const &c, room_t const &room, bool dim, bool dir) const;
	bool place_obj_along_wall(room_object type, room_t const &room, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
		cube_t const &place_area, unsigned objs_start, float front_clearance=0.0, bool add_door_clearance=0, unsigned pref_orient=4, bool pref_centered=0,
		colorRGBA const &color=WHITE, bool not_at_window=0, room_obj_shape shape=SHAPE_CUBE, float side_clearance=0.0, unsigned extra_flags=0, bool not_ext_wall=0, bool force_pref=0);
	bool place_model_along_wall(unsigned model_id, room_object type, room_t const &room, float height, rand_gen_t &rgen, float zval, unsigned room_id,
		float tot_light_amt, cube_t const &place_area, unsigned objs_start, float front_clearance=0.0, unsigned pref_orient=4,
		bool pref_centered=0, colorRGBA const &color=WHITE, bool not_at_window=0, unsigned extra_flags=0, bool force_pref=0, bool sideways=0);
	int check_valid_picture_placement(room_t const &room, cube_t const &c, float width, float zval, bool dim, bool dir, unsigned objs_start) const;
	void update_player_interact_objects(point const &player_pos);
	void update_creepy_sounds(point const &player_pos) const;
	point choose_creepy_sound_pos(point const &player_pos, rand_gen_t &rgen) const;
	void register_spark_floor_hit(point const &pos);
	bool line_intersect_walls(point const &p1, point const &p2, bool same_room=0) const;
	bool is_obj_pos_valid(room_object_t const &obj, bool keep_in_room, bool allow_block_door, bool check_stairs) const;
	bool is_rot_cube_visible(cube_t const &c, vector3d const &xlate, bool inc_mirror_reflections=0) const;
	bool is_cube_face_visible_from_pt(cube_t const &c, point const &p, unsigned dim, bool dir, bool same_room) const;
	bool check_obj_occluded(cube_t const &c, point const &viewer, occlusion_checker_noncity_t const &oc,
		bool reflection_pass=0, bool c_is_building_part=0, bool skip_basement_check=0) const;
	bool check_pg_br_wall_occlusion(point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area, vector3d const &view_dir) const;
	bool check_shelfrack_occlusion (point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area) const;
	bool check_warehouse_shelf_occlusion(point const &viewer, point const *const pts, unsigned npts, cube_t const &occ_area) const;
	bool is_entire_building_occluded(point const &viewer, occlusion_checker_noncity_t const &oc) const;
	bool register_indir_lighting_state_change(unsigned light_ix, bool is_door_change=0) const;
	bool is_attic_roof(tquad_with_ix_t const &tq, bool type_roof_only) const;
	bool has_ext_basement() const {return (interior && !interior->basement_ext_bcube.is_all_zeros());}
	bool has_pool        () const {return (interior &&  interior->pool.valid);}
	bool point_in_extended_basement(point const &pos) const;
	bool point_in_extended_basement_not_basement(point const &pos) const {return (point_in_extended_basement(pos) && !get_basement().contains_pt(pos));}
	bool cube_int_ext_basement(cube_t const &c) const {return (interior && interior->basement_ext_bcube.intersects(c));} // includes pools
	bool point_in_building_or_basement_bcube(point const &pos) const {return (bcube.contains_pt(pos) || point_in_extended_basement(pos));}
	bool point_in_extb_conn_room(point const &pos_bs) const;
	bool cube_intersects_extb_room(cube_t const &c, bool check_tunnel_pipes=0) const;
	bool cube_intersects_basement_or_extb_room(cube_t const &c, bool check_tunnel_pipes=0) const;
	bool point_in_courtyard(point const &pos_bs) const;
	bool point_on_basement_stairs(point const &pos_bs) const;
	bool is_cube_visible_through_extb_door(point const &viewer, cube_t const &c) const;
	float get_bcube_z1_inc_ext_basement() const {return (has_ext_basement() ? min(bcube.z1(), interior->basement_ext_bcube.z1()) : bcube.z1());}
	unsigned get_ext_basement_floor_ix(float zval) const;
	void get_pgbr_wall_ix_for_pos(point const &pos, index_pair_t &start, index_pair_t &end) const;
	cube_t get_bcube_inc_extensions () const;
	cube_t get_full_basement_bcube  () const;
	cube_t get_ext_basement_entrance() const;
	cube_t get_best_occluder(point const &camera_bs) const;
	cube_t get_step_for_ext_door(tquad_with_ix_t const &door) const;
	cube_t const &get_industrial_area() const {assert(is_industrial()); assert(!parts.empty()); return parts.front();}
	bool interior_visible_from_other_building_ext_basement(vector3d const &xlate, bool expand_for_light=0) const;
	bool top_of_mall_elevator_visible(point const &camera_bs, vector3d const &xlate) const;
	bool can_have_basement() const;
	void try_connect_ext_basement_to_building(building_t &b);
	void finalize_extb_conn_rooms(unsigned ds_start);
	template<typename T> void add_door_verts(cube_t const &D, T &drawer, door_rotation_t &drot, uint8_t door_type, bool dim, bool dir,
		float open_amt, bool opens_out, bool exterior, bool on_stairs=0, bool hinge_side=0, bool open_min_amt=0, bool draw_top_edge=0) const;
	tquad_with_ix_t set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj, bool exterior,
		float open_amt, bool opens_out, bool opens_up, bool swap_sides, bool open_min_amt, door_rotation_t &drot) const;
	tquad_with_ix_t set_interior_door_from_cube(door_t const &door) const;
	cube_t get_door_bounding_cube(door_t const &door) const;
	cube_t get_attic_access_door_avoid() const;
	cube_t get_light_switch_bounds(float floor_zval, float wall_edge, float wall_pos, bool dim, bool dir) const;
	void get_all_door_centers_for_room(room_t const &room, unsigned room_id, float zval, vector<point> &door_centers) const;
	void get_attic_windows(vect_tquad_with_ix_t &tquads, float offset_scale=1.0) const;
	void invalidate_nav_graph();
	void invalidate_nav_grid (unsigned floor_ix);
	point local_to_camera_space(point const &pos) const;
	void play_door_open_close_sound(point const &pos, bool open, float gain=1.0, float pitch=1.0) const;
	void play_open_close_sound(room_object_t const &obj, point const &sound_origin) const;
	void maybe_gen_chimney_smoke() const;
	int get_part_ix_containing_cube(cube_t const &c) const;
	int get_part_ix_containing_pt(point const &pt) const;
	cube_t get_part_containing_cube(cube_t const &c) const;
	cube_t get_part_containing_pt(point const &pt) const;
	bool move_sphere_to_valid_part(point &pos, point const &p_last, float radius) const;
	void remove_paint_in_cube(cube_t const &c) const;
	bool has_water() const {return (interior && interior->water_zval != 0.0);}
	bool water_visible_to_player() const;
	float get_floor_below_water_level() const;
	cube_t get_water_cube(bool full_room_height=0) const;
	bool point_in_water_area(point const &p, bool full_room_height=1) const;
	bool point_in_or_above_pool(point const &pt) const;
	bool set_float_height(point &pos, float radius, float ceil_zval, float density=0.5) const;
	bool find_mirror_needing_reflection(vector3d const &xlate) const;
	float get_elevator_floor_spacing(elevator_t            const &e) const {return ( e.in_mall       ? get_mall_floor_spacing() : get_window_vspace());}
	float get_stairs_floor_spacing  (stairs_landing_base_t const &s) const {return ((s.in_mall == 1) ? get_mall_floor_spacing() : get_window_vspace());}
	float get_room_floor_spacing    (room_t                const &r) const {return ( r.is_mall()     ? get_mall_floor_spacing() : get_window_vspace());}
	void print_building_manifest() const;
	void print_building_stats() const;
private:
	void create_per_part_ext_verts();
	void finish_gen_geometry(rand_gen_t &rgen, bool has_overlapping_cubes);
	void assign_name(rand_gen_t &rgen);
	std::string get_name_for_floor(unsigned floor_ix) const;
	bool add_outdoor_ac_unit(rand_gen_t &rgen);
	bool add_chimney(bool two_parts, bool stacked_parts, bool hipped_roof[4], float roof_dz[4], unsigned force_dim[2], rand_gen_t &rgen);
	float get_min_hallway_width() const;
	bool can_use_hallway_for_part(unsigned part_id) const;
	cube_t get_hallway_for_part(cube_t const &part, float &num_hall_windows, float &hall_width, float &room_width);
	void gen_interior_int(rand_gen_t &rgen, bool has_overlapping_cubes);
	void maybe_add_basement(rand_gen_t rgen);
	bool extend_underground_basement(rand_gen_t rgen);
	bool is_basement_room_not_int_bldg(cube_t const &room, building_t const *exclude=nullptr, bool allow_outside_grid=0) const;
	bool is_basement_room_under_mesh_not_int_bldg(cube_t const &room, building_t const *exclude=nullptr, bool allow_outside_grid=0) const;
	bool is_basement_room_placement_valid(cube_t &room, ext_basement_room_params_t &P, bool dim, bool dir, bool *add_end_door=nullptr, building_t const *exclude=nullptr) const;
	bool add_underground_exterior_rooms(rand_gen_t &rgen, cube_t const &door_bcube, cube_t const &basement, bool wall_dim, bool wall_dir, float length_mult);
	bool check_pool_room_slice_valid(cube_t const &slice, int skip_room_ix) const;
	void maybe_assign_extb_room_as_swimming(rand_gen_t &rgen);
	void add_wall_section_above_pool_room_door(door_stack_t &ds, room_t const &room);
	unsigned setup_multi_floor_room(extb_room_t &room, door_t const &door, bool wall_dim, bool wall_dir, rand_gen_t &rgen);
	void add_backrooms_droplet_spawners(rand_gen_t rgen);
	void update_droplet_spawners();
	bool add_ext_basement_rooms_recur(extb_room_t &parent_room, ext_basement_room_params_t &P, float door_width, bool dim, unsigned depth, rand_gen_t &rgen);
	unsigned max_expand_underground_room(cube_t &room, bool dim, bool dir, bool is_mall, rand_gen_t &rgen) const;
	void add_mall_se_landing(cube_t const &c, bool is_escalator, bool se_dim, bool se_dir, bool ww_dir);
	void setup_mall_concourse(cube_t const &room, bool dim, bool dir, rand_gen_t &rgen);
	bool is_cube_city_placement_invalid(cube_t const &c) const;
	bool is_store_placement_invalid(cube_t const &store) const;
	void add_mall_stores(cube_t const &room, bool dim, bool entrance_dir, rand_gen_t &rgen);
	void add_mall_store(cube_t const &store, cube_t const &window_area, bool dim, bool dir, bool &has_adj_store, rand_gen_t &rgen);
	void add_extb_room_floor_and_ceil(cube_t const &room);
	void add_mall_stairs();
	bool adjust_zval_for_mall_stairs(point const &pos, float &zval) const;
	float get_mall_floor_spacing(cube_t const &room) const;
	bool inside_mall_hallway(point const &pos) const;
	bool is_inside_mall_stores(point const &pos) const;
	room_t const &get_mall_concourse() const {assert(has_mall()); return interior->get_extb_start_room();}
	float get_mall_floor_spacing    () const {return get_mall_floor_spacing(get_mall_concourse());}
	static float get_mall_top_window_gap(float mall_floor_spacing, float window_vspace) {return 0.5*(mall_floor_spacing - window_vspace);}
	cube_t get_mall_center(cube_t const &room) const;
	void get_mall_open_areas(cube_t const &room, vect_cube_t &openings) const;
	cube_t add_ext_basement_door(cube_t const &room, float door_width, bool dim, bool dir, bool is_end_room, bool is_tall_room, rand_gen_t &rgen, bool opens_other_side=0);
	cube_t add_and_connect_ext_basement_room(extb_room_t &room, ext_basement_room_params_t &P,
		float door_width, bool dim, bool dir, bool is_end_room, unsigned depth, bool const add_doors[2], rand_gen_t &rgen);
	void end_ext_basement_hallway(extb_room_t &room, cube_t const &conn_bcube, ext_basement_room_params_t &P,
		float door_width, bool dim, bool dir, unsigned depth, rand_gen_t &rgen);
	void get_valid_extb_room_end_doors(room_t const &room, float zval, unsigned room_id, float end_pad_ext, cube_with_ix_t doors[2]) const;
	void add_false_door_to_extb_room_if_needed(room_t const &room, float zval, unsigned room_id);
	bool is_tunnel_bcube_placement_valid(cube_t const &tunnel_bc) const;
	bool is_tunnel_placement_valid(point const &p1, point const &p2, float radius) const;
	void try_extend_tunnel(point &p1, point &p2, float max_extend, float check_radius, unsigned dim, bool dir) const;
	bool try_place_tunnel_at_extb_hallway_end(room_t &room, unsigned room_id, rand_gen_t &rgen);
	bool add_extend_tunnel_seg(point const &end_pt, float radius, float wall_thickness, float min_len, float max_extend, float gate_dist_from_end,
		unsigned parent_seg_ix, unsigned seg_depth, bool dim, bool pdir, bool cdir, rand_gen_t &rgen);
	building_t *get_conn_bldg_for_pt(point const &p, float radius=0.0) const;
	building_t *get_bldg_containing_pt(point const &p);
	bool is_visible_through_conn(building_t const &b, vector3d const &xlate, float view_dist, bool expand_for_light=0) const;
	cube_t get_conn_room_closest_to(point const &pos_bs) const;
	bool has_L_shaped_roof_area() const;
	void get_attic_roof_tquads(vect_tquad_with_ix_t &tquads) const;
	bool add_attic_access_door(cube_t const &ceiling, unsigned part_ix, unsigned num_floors, unsigned rooms_start, rand_gen_t &rgen);
	bool is_light_placement_valid(cube_t const &light, room_t const &room, float pad) const;
	void try_place_light_on_ceiling(cube_t const &light, room_t const &room, unsigned room_id, bool room_dim, float pad, bool allow_rot,
		bool allow_mult, unsigned nx, unsigned ny, unsigned check_coll_start, vect_cube_t &lights, rand_gen_t &rgen) const;
	void try_place_light_on_wall   (cube_t const &light, room_t const &room, bool room_dim, float zval, vect_cube_t &lights, rand_gen_t &rgen) const;
	bool clip_cube_to_parts(cube_t &c, vect_cube_t &cubes) const;
	cube_t get_walkable_room_bounds(room_t const &room, bool floor_space_only=0) const;
	bool is_cube_contained_in_parts(cube_t const &c) const;
	void expand_ground_floor_cube(cube_t &cube, cube_t const &skip=cube_t()) const;
	void get_exclude_cube(point const &pos, cube_t const &skip, cube_t &exclude, bool camera_in_building) const;
	void move_door_to_other_side_of_wall(tquad_with_ix_t &door, float dist_mult, bool invert_normal) const;
	void clip_door_to_interior(tquad_with_ix_t &door) const;
	void cut_holes_for_ext_doors(building_draw_t &bdraw, point const &contain_pt, unsigned draw_parts_mask, cube_t const &clamp_cube=cube_t()) const;
	bool is_valid_door_pos(cube_t const &door, float door_width, bool dim) const;
	bool is_cube_close_to_exterior_doorway(cube_t const &c, float dmin=0.0, bool inc_open=0) const;
	bool is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin=0.0, bool inc_open=0, bool check_open_dir=0) const;
	bool is_obj_placement_blocked(cube_t const &c, cube_t const &room, bool inc_open_doors, bool check_open_dir=0, float dmin=0.0) const;
	bool is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, vector2d const &room_pad=vector2d()) const;
	bool check_cube_intersect_walls(cube_t const &c) const;
	bool check_cube_contained_in_part(cube_t const &c) const;
	bool is_valid_stairs_elevator_placement(cube_t const &c, cube_t const &c_nopad, float pad, int dim=2, bool check_walls=1, bool check_private_rooms=0) const;
	bool clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const;
	void add_ceiling_cube_no_skylights(cube_t const &c);
	void calc_room_ext_sides(room_t &room) const;
	unsigned add_room(cube_t const &room, unsigned part_id, unsigned num_lights=1, bool is_hallway=0, bool is_office=0, bool is_sec_bldg=0);
	void add_or_extend_elevator(elevator_t const &elevator, bool add);
	void remove_intersecting_roof_cubes(cube_t const &c);
	bool overlaps_other_room_obj(cube_t const &c, unsigned objs_start=0, bool check_all=0, unsigned const *objs_end=nullptr) const;
	bool overlaps_obj_or_placement_blocked(cube_t const &c, cube_t const &room, unsigned objs_start, bool check_all=0, float dmin=0.0) const;
	bool overlaps_any_placed_obj(cube_t const &c) const;
	bool overlaps_or_adj_int_window (cube_t const &c) const;
	bool check_skylight_intersection(cube_t const &c) const;
	int classify_room_wall(room_t const &room, float zval, bool dim, bool dir, bool ret_sep_if_part_int_part_ext) const;
	unsigned count_ext_walls_for_room(room_t const &room, float zval) const;
	bool room_has_stairs_or_elevator(room_t const &room, float zval, unsigned floor) const;
	bool is_room_office_bathroom(room_t &room, float zval, unsigned floor) const;
	int gather_room_placement_blockers(cube_t const &room, unsigned objs_start, vect_cube_t &blockers, bool inc_open_doors=1, bool ignore_chairs=0) const;
	void get_doorways_for_room(cube_t const &room, float zval, vect_door_stack_t &doorways, bool all_floors=0) const;
	vect_door_stack_t &get_doorways_for_room(cube_t const &room, float zval, bool all_floors=0) const;
	bool is_room_an_exit(cube_t const &room, int room_ix, float zval) const;
	vector3d get_office_chair_size() const;
	bool add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos, colorRGBA const &chair_color, bool dim,
		bool dir, float tot_light_amt, bool office_chair=0, bool enable_rotation=0, bool bar_stool=0, bool no_push_out=0, int wooden_or_plastic=2, bool reduced_clearance=0);
	unsigned add_table_and_chairs(rand_gen_t rgen, room_t const &room, vect_cube_t &blockers, unsigned room_id, point const &place_pos,
		colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, unsigned max_chairs=4, bool use_tall_table=0, int wooden_or_plastic=2);
	void shorten_chairs_in_region(cube_t const &region, unsigned objs_start);
	void add_trashcan_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool check_last_obj);
	bool add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	bool add_desk_to_room    (rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, colorRGBA const &chair_color, float zval,
		unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement, unsigned desk_ix=0, bool no_computer=0, bool force_computer=0, bool add_phone=0);
	void add_papers_to_surface(cube_t const &c, bool dim, bool dir, unsigned max_num, rand_gen_t &rgen, unsigned room_id, float tot_light_amt, cube_t const &avoid=cube_t());
	void add_pens_pencils_to_surface(cube_t const &c, bool dim, bool dir, unsigned max_num, rand_gen_t &rgen, unsigned room_id, float tot_light_amt, cube_t const &avoid=cube_t());
	void add_filing_cabinet_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_office_objs     (rand_gen_t rgen, room_t const &room, vect_cube_t &blockers, colorRGBA const &chair_color,
		float zval, unsigned room_id, unsigned floor, float tot_light_amt, unsigned objs_start, bool is_basement);
	void add_clock_to_room_wall(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, int type=2);
	void add_industrial_office_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool create_office_cubicles(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	void add_office_pillars(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t const &lights, vect_cube_t &blockers);
	vector2d get_conf_room_table_length_width(cube_t const &room) const;
	bool add_conference_objs (rand_gen_t rgen,  room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned floor_ix);
	void add_lounge_objs     (rand_gen_t rgen,  room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_lobby);
	bool add_vending_machine (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t const &place_area);
	bool add_wall_tv(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool check_valid_closet_placement(cube_t const &c, room_t const &room, unsigned objs_start, unsigned bed_ix, float min_bed_space=0.0) const;
	bool add_bedroom_objs    (rand_gen_t rgen, room_t &room, vect_cube_t &blockers, colorRGBA const &chair_color, float zval, unsigned room_id, unsigned floor,
		float tot_light_amt, unsigned objs_start, bool room_is_lit, bool is_basement, bool force, light_ix_assign_t &light_ix_assign);
	bool replace_light_with_ceiling_fan(rand_gen_t &rgen, cube_t const &room, cube_t const &avoid, unsigned room_id, float tot_light_amt, unsigned light_obj_ix);
	bool add_bed_to_room     (rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval,
		unsigned room_id, float tot_light_amt, unsigned floor, bool force, int &bed_size_ix, room_object_t const &other_bed);
	bool add_ball_to_room    (rand_gen_t &rgen, room_t const &room, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
		int force_type=-1, cube_t const &avoid_xy=cube_t(), bool in_pool=0);
	bool maybe_add_fireplace_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t &blockers, float zval, unsigned room_id, float tot_light_amt);
	bool add_classroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start,
		colorRGBA const &chair_color, unsigned &td_orient);
	bool add_classroom_desk(rand_gen_t &rgen, room_t const &room, cube_t const &desk, unsigned room_id, float tot_light_amt,
		colorRGBA const &chair_color, bool dim, bool dir, unsigned desk_ix);
	void add_objects_next_to_classroom_chalkboard(rand_gen_t &rgen, room_object_t const &cb, room_t const &room, float zval, unsigned objs_start);
	void add_hallway_lockers(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool get_hospital_room_bathroom(room_t const &room, unsigned room_id, int &nested_room_ix, cube_t &bathroom) const;
	bool add_hospital_room_objs (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start, int &nested_room_ix);
	bool add_waiting_room_objs  (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start, int &nested_room_ix);
	void place_chairs_along_walls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
		colorRGBA const &chair_color, bool is_plastic, unsigned num);
	bool add_exam_room_objs     (rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start);
	bool add_operating_room_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, float &tot_light_amt, unsigned objs_start, unsigned lights_start);
	void add_hospital_medicine_cabinet(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	float add_flooring       (room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned flooring_type);
	bool add_bathroom_objs   (rand_gen_t rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt,
		unsigned objs_start, unsigned floor, bool is_basement, bool add_shower_tub, unsigned &added_bathroom_objs_mask);
	bool add_vanity_to_room  (rand_gen_t &rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_bathroom_plumbing(room_object_t const &obj);
	bool add_tp_roll         (cube_t const &room, unsigned room_id, float tot_light_amt, bool dim, bool dir, float length, float zval, float wall_pos, bool check_valid_pos=0);
	bool divide_bathroom_into_stalls(rand_gen_t &rgen, room_t &room, float zval, unsigned room_id, float tot_light_amt, unsigned floor);
	void add_door_sign                (std::string const &text, room_t const &room, float zval, unsigned room_id, bool no_check_adj_walls=0, cube_t const &avoid=cube_t());
	void add_numbered_door_sign       (std::string const &text, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, bool no_check_adj_walls=0);
	void add_door_sign_remove_existing(std::string const &text, room_t const &room, float zval, unsigned room_id, unsigned objs_start);
	void add_office_door_sign(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id);
	void add_out_or_order_sign(cube_t const &door_bc, bool dim, bool dir, unsigned room_id);
	void make_door_out_or_order(room_t const &room, float zval, unsigned room_id, unsigned door_stack_ix);
	bool add_kitchen_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool allow_adj_ext_door);
	bool add_fishtank_to_room(rand_gen_t&rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t const &place_area);
	bool add_livingroom_objs (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_diningroom_objs (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_library_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	bool add_storage_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement, bool has_stairs);
	void add_boxes_and_crates(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
		unsigned num_max, bool is_basement, cube_t const &room_bounds, cube_t const &crate_bounds, vect_cube_t const &exclude);
	bool add_interrogation_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_buckets_to_room (rand_gen_t &rgen, cube_t place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num);
	bool add_ladder_to_room  (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_garage_objs     (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	void add_floor_clutter_objs(rand_gen_t  rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_floor_clutter_objs(rand_gen_t &rgen, room_t const &room, cube_t place_area, float zval, unsigned room_id,
		float tot_light_amt, unsigned objs_start, bool add_bottles, bool add_trash, bool add_papers, bool add_glass);
	void add_basement_clutter_objs(rand_gen_t  rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	unsigned add_water_heaters (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool single_only=0);
	bool add_basement_utility_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_furnace_to_room (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_breaker_panel   (rand_gen_t &rgen, cube_t const &c, float ceil_zval, bool dim, bool dir, unsigned room_id, float tot_light_amt);
	bool add_office_utility_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_server_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_security_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_laundry_basket  (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t place_area);
	bool add_laundry_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned &added_bathroom_objs_mask);
	void add_couches_to_room (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned const counts[4]);
	bool add_pool_room_objs  (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	void add_swimming_pool_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	void add_retail_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, light_ix_assign_t &light_ix_assign);
	float gather_room_lights(unsigned objs_start_inc_lights, vect_cube_t &lights) const;
	void add_ladders_to_nested_room_roofs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float light_amt, cube_t const &place_area);
	cube_t add_factory_ladders_and_catwalks(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float light_amt, cube_t const &place_area,
		cube_t const &place_area_upper, float beams_z1, float support_width, vect_cube_t supports[2], vect_cube_t const &lights, vect_cube_t &ladders, vector<float> const &beam_pos);
	void add_industrial_ducts_and_hvac(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float light_amt, float support_width, float ceil_zval,
		cube_t const &place_area_upper, vect_cube_t const &beams, vect_cube_t supports[2], unsigned objs_start);
	void add_industrial_sprinkler_pipes(rand_gen_t &rgen, room_t const &room, unsigned room_id,
		float support_width, unsigned objs_start, vect_cube_t const &lights, vect_cube_t const &beams);
	void add_warehouse_shelves(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start,
		float light_amt, float support_width, cube_t const &place_area, vect_cube_t const &supports);
	void add_industrial_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start_inc_lights);
	void add_retail_pillar(cube_t const &pillar, float zval, unsigned room_id, bool is_tall);
	void add_U_stair_landing_lights(stairwell_t const &s, unsigned room_id, unsigned light_ix, float floor_zval);
	void add_checkout_objs   (cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool dim, bool dir, bool cr_dir);
	void add_shelves(cube_t const &c, bool dim, bool dir, unsigned room_id, float tot_light_amt, unsigned flags, unsigned item_flags, rand_gen_t &rgen);
	void add_shelf_rack(cube_t const &c, bool dim, unsigned style_id, unsigned &rack_id, unsigned room_id,
		unsigned extra_flags, unsigned item_category, bool add_occluders, rand_gen_t &rgen);
	bool maybe_add_walkway_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, light_ix_assign_t &light_ix_assign);
	void add_clock(cube_t const &clock, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool digital);
	void add_clock_to_cube(cube_t const &c, float zval, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool digital);
	bool add_fire_ext_along_wall(cube_t const &room, float zval, unsigned room_id, float tot_light_amt, bool dim, bool dir, rand_gen_t &rgen);
	void add_fire_ext        (float height, float radius, float zval, float wall_edge, float pos_along_wall,
		unsigned room_id, float tot_light_amt, bool dim, bool dir, bool center_mount=0);
	bool is_contained_in_wall_range(float wall_pos, float cov_lo, float cov_hi, float zval, bool dim) const;
	void add_pri_hall_objs   (rand_gen_t rgen, rand_gen_t room_rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned floor_ix, unsigned objs_start);
	void add_wall_us_flag    (float wall_pos, float flag_pos, float zval, bool dim, bool dir, unsigned room_id, float tot_light_amt);
	bool add_reception_desk  (rand_gen_t &rgen, cube_t const &desk, bool dim, bool dir, unsigned room_id, float tot_light_amt);
	void add_cameras_to_room (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void assign_attic_type   (rand_gen_t rgen);
	void add_attic_objects   (rand_gen_t rgen);
	void add_attic_ductwork  (rand_gen_t rgen, room_object_t const &furnace, vect_cube_t &avoid_cubes);
	bool add_attic_roof_vent(point const &bot_center, float radius, unsigned room_id, float light_amt=1.0);
	int choose_air_intake_room() const;
	int vent_in_attic_test(cube_t const &vent, bool dim) const;
	void add_exterior_ac_pipes(rand_gen_t rgen);
	void add_tunnel_objects   (rand_gen_t rgen);
	void add_interior_window_objects();
	void add_padlocks(rand_gen_t rgen);
	bool add_padlock_to_door     (unsigned door_ix, unsigned lock_color_mask, rand_gen_t &rgen);
	bool remove_padlock_from_door(unsigned door_ix, point const &remove_pos);
	vector3d get_parked_car_size() const;
	cube_t get_ext_basement_door_blocker() const;
	void add_parking_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
		unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z, light_ix_assign_t &light_ix_assign);
	void add_backrooms_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light);
	unsigned add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, vect_cube_t &rooms_to_light);
	void add_tv_to_wall(cube_t const &tv, unsigned room_id, float light_amt, bool dim, bool dir, bool use_monitor_image=0, int on_off=2);
	bool is_valid_placement_at_mall_wall(cube_t const &c, cube_t const &entrance_stairs_bcube, vect_cube_t const &wall_blockers) const;
	void add_large_trashcan_contents(rand_gen_t &rgen, room_object_t const &tcan, unsigned room_id, float tot_light_amt);
	bool add_mall_table_with_chairs(rand_gen_t &rgen, cube_t const &table, cube_t const &place_area, colorRGBA const &chair_color,
		unsigned room_id, float tot_light_amt, bool dim, unsigned tid_tag, vect_cube_t &blockers);
	bool add_food_court_objs(rand_gen_t &rgen, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, vect_cube_t const &blockers);
	void add_mall_store_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned &type_mask, light_ix_assign_t &light_ix_assign);
	void add_ceiling_ducts(cube_t const &room, float ceil_zval, unsigned room_id, bool dim, unsigned skip_dir, float light_amt,
		bool cylin_ducts, bool skip_ends, bool skip_top, rand_gen_t &rgen);
	void add_row_of_bookcases(cube_t const &row, float zval, unsigned room_id, float light_amt, bool dim, bool place_inside);
	void add_shelves_along_walls(cube_t const &room_area, float zval, unsigned room_id, float light_amt,
		bool dim, unsigned store_type, float height, float depth, bool place_inside, rand_gen_t &rgen);
	void add_missing_backrooms_lights(rand_gen_t rgen, float zval, unsigned room_id, unsigned objs_start, unsigned lights_start,
		room_object_t const &ref_light, vect_cube_t const &rooms_to_light, light_ix_assign_t &light_ix_assign);
	void add_mall_lower_floor_lights(room_t const &room, unsigned room_id, unsigned lights_start, light_ix_assign_t &light_ix_assign);
	void add_sub_room_light(room_object_t light, room_t const &room, unsigned room_id, bool dim, unsigned objs_start, light_ix_assign_t &light_ix_assign, rand_gen_t &rgen);
	bool add_basement_pipes(vect_cube_t const &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, vect_riser_pos_t const &risers, vect_cube_t &pipe_cubes,
		unsigned room_id, unsigned num_floors, unsigned objs_start, float ceil_zval, rand_gen_t &rgen, unsigned pipe_type, bool allow_place_fail=0);
	void add_ext_basement_hallway_pipes_recur(unsigned room_id, bool hall_dim, unsigned pipe_type, float radius_factor,
		pipe_t const &parent, vector<pipe_t> &pipes, vector<pipe_t> &fittings, rand_gen_t &rgen) const;
	bool add_sprinkler_pipes(vect_cube_t const &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, vect_cube_t const &pipe_cubes,
		unsigned room_id, unsigned num_floors, unsigned objs_start, rand_gen_t &rgen,
		float custom_floor_spacing=0.0, float wall_pad=0.0, unsigned pref_dim=2, vect_cube_t const &vpipe_avoid=vect_cube_t());
	void get_pipe_basement_water_connections(vect_riser_pos_t &sewer, vect_riser_pos_t &cold_water, vect_riser_pos_t &hot_water, rand_gen_t &rgen) const;
	void get_pipe_basement_gas_connections  (vect_riser_pos_t &pipes) const;
	void add_basement_electrical(vect_cube_t &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, int room_id, rand_gen_t &rgen);
	void add_basement_electrical_house(rand_gen_t &rgen);
	void add_house_basement_pipes(rand_gen_t &rgen);
	void place_book_on_obj   (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, unsigned objs_start, bool use_dim_dir);
	bool place_bottle_on_obj (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t(), bool place_at_z1=0);
	bool place_dcan_on_obj   (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t(), bool place_at_z1=0);
	bool place_plant_on_obj  (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, float sz_scale, vect_cube_t const &avoid=vect_cube_t());
	bool place_laptop_on_obj (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t(), bool use_dim_dir=0);
	bool place_pizza_on_obj  (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t());
	bool place_plate_on_obj  (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t());
	bool place_cup_on_obj    (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t(), bool make_empty=0);
	bool place_toy_on_obj    (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t());
	bool place_banana_on_obj (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid=vect_cube_t());
	bool place_phone_on_obj  (rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, bool dim, bool dir, float overhang_amt=0.0);
	bool add_rug_to_room     (rand_gen_t  rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool hang_pictures_whiteboard_chalkboard_in_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
		unsigned floor_ix, bool is_basement, unsigned pref_orient=4);
	void add_plants_to_room  (rand_gen_t  rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num);
	void add_boxes_to_room   (rand_gen_t  rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned max_num);
	void add_stains_to_room  (rand_gen_t  rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_floor_stains    (rand_gen_t &rgen, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num, float rmax);
	void add_light_switches_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor, bool is_basement);
	void add_outlets_to_room  (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor, bool is_basement, bool is_kitchen=0);
	bool add_wall_vent_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool check_for_ducts);
	bool add_ceil_vent_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start);
	bool check_if_placed_on_interior_wall(cube_t const &c, room_t const &room, bool dim, bool dir) const;
	bool place_eating_items_on_table(rand_gen_t &rgen, unsigned table_obj_id);
	void place_objects_onto_surfaces(rand_gen_t rgen, room_t const &room, unsigned room_id, float tot_light_amt,
		unsigned objs_start, unsigned floor, bool is_basement, bool not_private);
	bool can_be_bedroom_or_bathroom(room_t const &room, unsigned floor_ix, bool skip_conn_check=0) const;
	bool can_be_bathroom(room_t const &room) const;
	bool find_mirror_in_room(unsigned room_id, vector3d const &xlate, float &dmin_sq, bool same_room) const;
	int find_main_roof_tquad_ix(rand_gen_t &rgen, bool skip_if_has_other_obj) const;
	void add_chimney_cap(rand_gen_t &rgen);
	void maybe_add_fire_escape(rand_gen_t &rgen);
	void add_balconies(rand_gen_t &rgen, vect_cube_t &balconies);
	void add_gutter_downspouts(rand_gen_t &rgen, vect_cube_t const &balconies);
	void add_extra_obj_slots();
	void add_wall_and_door_trim_if_needed();
	void add_trim_for_door_or_int_window(cube_t const &c, bool dim, bool draw_top_edge, bool draw_bot_trim,
		float side_twidth, float top_twidth, float side_texp, float floor_spacing, float extra_top_gap=0.0);
	void add_wall_and_door_trim();
	void add_window_trim_and_coverings(bool add_trim, bool add_blinds, bool add_ext_sills=0);
	void add_pool_trim();
	bool check_ext_step_valid(cube_t const &c, unsigned ext_objs_start, unsigned exclude_ix, float head_clearance) const;
	void add_ext_door_steps(unsigned ext_objs_start);
	unsigned count_num_int_doors(room_t const &room) const;
	bool check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs) const;
	void split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen);
	bool test_coll_with_sides(point &pos, point const &p_last, float radius, vector3d const &xlate, cube_t const &part, unsigned part_id, vector3d *cnorm) const;
	void gather_interior_cubes(vect_colored_cube_t &cc, cube_t const &ext_bcube) const;
	void get_lights_with_priorities(point const &target, cube_t const &valid_area, vector<pair<float, unsigned>> &lights_to_sort) const;
	void get_all_windows(vect_cube_with_ix_t &windows) const;
	void register_indir_lighting_geom_change() const;
	void register_blinds_state_change() const;
	bool is_light_occluded(point const &lpos, point const &camera_bs) const;
	void clip_ray_to_walls(point const &p1, point &p2, vect_cube_t const walls[2]) const;
	void refine_light_bcube(point const &lpos, float light_radius, room_t const &room, cube_t &light_bcube, bool is_parking_garage) const;
	cube_t get_rotated_bcube(cube_t const &c, bool inv_rotate=0) const;
	cube_t const &get_part_for_room(room_t const &room) const {assert(room.part_id < parts.size()); return parts[room.part_id];}
	bool are_parts_stacked(cube_t const &p1, cube_t const &p2) const;
	room_type get_room_type_and_floor(int room_id, float zval, unsigned &floor_ix) const;
	void add_window_coverings(cube_t const &window, bool dim, bool dir);
	void add_window_blinds(cube_t const &window, bool dim, bool dir, unsigned room_ix, unsigned floor);
	void add_bathroom_window(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor);
	int get_room_id_for_window(cube_t const &window, bool dim, bool dir, bool &is_split) const;
	void register_open_ext_door_state(int door_ix);
	void add_interior_door(door_t &door, bool is_bathroom=0, bool make_unlocked=0, bool make_closed=0);
	void add_interior_door_for_floor(door_t &door, bool is_bathroom=0, bool make_unlocked=0, bool make_closed=0);
	void remove_section_from_cube_and_add_door(cube_t &c, cube_t &c2, float v1, float v2, bool xy,
		bool open_dir, bool is_bathroom=0, bool make_unlocked=0, bool make_closed=0);
	void insert_door_in_wall_and_add_seg(cube_t &wall, float v1, float v2, bool dim, bool open_dir,
		bool keep_high_side=0, bool is_bathroom=0, bool make_unlocked=0, bool make_closed=0);
	void reverse_door_hinges_if_needed();
	void ensure_doors_to_room_are_closed(room_t const &room, unsigned doors_start, bool ensure_locked=0);
	unsigned get_floor_for_zval(float zval) const {return unsigned((zval - get_bcube_z1_inc_ext_basement())/get_window_vspace());}
	building_loc_t get_building_loc_for_pt(point const &pos) const;
	bool same_room_and_floor_as_player(person_t const &person) const;
	bool is_player_visible(person_t const &person, unsigned vis_test) const;
	bool can_target_player(person_t const &person) const;
	bool need_to_update_ai_path(person_t &person) const;
	void set_bcube_from_rotated_cube(cube_t const &bc);
	bool apply_paint(point const &pos, vector3d const &dir, colorRGBA const &color, unsigned emissive_color_id, room_object const obj_type) const;
	bool apply_toilet_paper(point const &pos, vector3d const &dir, float half_width);
	void register_button_event(room_object_t const &button);
	void call_elevator_to_floor(elevator_t &elevator, unsigned floor_ix, bool is_inside_elevator, bool is_up);
	void call_elevator_to_floor_and_light_nearest_button(elevator_t &elevator, unsigned floor_ix, bool is_inside_elevator, bool is_up);
	void run_ball_update(vect_room_object_t::iterator ball_it, point const &player_pos, float player_z1, bool player_is_moving);
	void update_pool_table(room_object_t &ball);
	bool get_zval_for_pool_bottom(point const &pos, float &zval) const;
	bool get_zval_of_floor(point const &pos, float radius, float &zval) const;
	bool get_zval_for_obj_placement(point const &pos, float radius, float &zval, bool add_z_bias) const;
	void register_player_death(point const &camera_bs);
	void add_blood_decal(point const &pos, float radius, colorRGBA const &color=WHITE);
	void add_broken_glass_to_floor(point const &pos, float radius);
	void add_broken_glass_decal(point const &pos, float radius, rand_gen_t &rgen);
	void play_tape_sound(point const &sound_pos, float sound_gain, bool tape_break) const;
	bool is_obj_above_ramp(cube_t const &c) const;
	bool is_room_above_ramp(cube_t const &room, float zval) const;
	bool is_room_lit(int room_id, float zval) const;
	void get_rooms_for_door(unsigned door_ix, int room_ix[2]) const;
	void get_lights_for_room_and_floor(unsigned room_ix, unsigned floor_ix, vector<unsigned> &light_ids) const;
	void get_lights_near_door(unsigned door_ix, vector<unsigned> &light_ids) const;
	bool is_cube_visible_through_door(point const &viewer, cube_t const &c, door_t const &door) const;
	void set_rgen_state_for_building(rand_gen_t &rgen) const;
public:
	// ray queries
	bool check_line_intersect_doors(point const &p1, point const &p2, bool inc_open=0) const;
	bool is_pt_visible(point const &p1, point const &p2) const;
	bool is_sphere_visible(point const &center, float radius, point const &pt) const;
	bool is_pt_lit(point const &pt) const;
	bool is_sphere_lit(point const &center, float radius) const;
}; // end building_t

struct vect_building_t : public vector<building_t> {
	void ai_room_update(float delta_dir, float dmax, point const &camera_bs, rand_gen_t &rgen);
};

struct building_draw_utils {
	static void calc_normals(building_geom_t const &bg, vector<vector3d> &nv, unsigned ndiv);
	static void calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, cube_t const &part, vect_point &pts);
};

struct walkway_base_t {
	bool open_ends[2]={};
	float floor_spacing;
	unsigned side_mat_ix, roof_mat_ix; // matches building material
	colorRGBA side_color, roof_color;
	walkway_base_t(unsigned smix, unsigned rmix, colorRGBA const &sc, colorRGBA const &rc, float floor_spacing_) :
		floor_spacing(floor_spacing_), side_mat_ix(smix), roof_mat_ix(rmix), side_color(sc), roof_color(rc) {}
};
struct bldg_walkway_t : public cube_t, public walkway_base_t {
	bool dim;
	bldg_walkway_t(cube_t const &c, bool d, unsigned smix, unsigned rmix, colorRGBA const &sc, colorRGBA const &rc, float fs) :
		cube_t(c), walkway_base_t(smix, rmix, sc, rc, fs), dim(d) {}
};
typedef vector<bldg_walkway_t> vect_bldg_walkway_t;

class city_lights_manager_t {
protected:
	cube_t lights_bcube;
	float light_radius_scale=1.0, dlight_add_thresh=0.0;
	bool prev_had_lights=0;
public:
	virtual ~city_lights_manager_t() {}
	cube_t get_lights_bcube() const {return lights_bcube;}
	void add_player_flashlight(float radius_scale);
	void tighten_light_bcube_bounds(vector<light_source> const &lights);
	void clamp_to_max_lights(vector3d const &xlate, vector<light_source> &lights);
	bool begin_lights_setup(vector3d const &xlate, float light_radius, vector<light_source> &lights);
	void finalize_lights(vector<light_source> &lights);
	void setup_shadow_maps(vector<light_source> &light_sources, point const &cpos, unsigned max_smaps, bool sec_camera_mode=0);
	virtual bool enable_lights() const = 0;
};

struct ped_draw_vars_t {
	building_t &building; // Note: building ref is non-const because we may update rendering data for the people inside it
	occlusion_checker_noncity_t &oc;
	shader_t &s;
	vector3d const &xlate;
	unsigned bix;
	bool shadow_only, reflection_pass, in_retail_room;

	ped_draw_vars_t(building_t &b, occlusion_checker_noncity_t &oc_, shader_t &s_, vector3d const &x, unsigned bix_, bool so, bool rp, bool irr=0)
		: building(b), oc(oc_), s(s_), xlate(x), bix(bix_), shadow_only(so), reflection_pass(rp), in_retail_room(irr) {}
};

class water_sound_manager_t {
	point const camera_bs;
	point closest; // in camera space
	float dmin_sq=0.0;
public:
	water_sound_manager_t(point const &camera_bs_) : camera_bs(camera_bs_) {}
	void register_running_water(room_object_t const &obj, building_t const &building);
	void finalize();
};

inline void clip_low_high_tc(float &t0, float &t1) {
	if (fabs(t0 - t1) < 0.5) {t0 = t1 = 0.0;} // too small to have a window
	else {t0 = round_fp(t0); t1 = round_fp(t1);} // Note: round() is much faster than nearbyint(), and round_fp() is faster than round()
}
inline bool check_bcube_sphere_coll(cube_t const &bcube, point const &sc, float radius, bool xy_only) {
	return (xy_only ? sphere_cube_intersect_xy(sc, radius, bcube) : sphere_cube_intersect(sc, radius, bcube));
}

template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &cubes, unsigned start_ix) { // T must derive from cube_t
	assert(start_ix <= cubes.size());
	for (auto i = cubes.begin()+start_ix; i < cubes.end(); ++i) {if (i->intersects(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_no_adj(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_xy(cube_t const &bcube, vector<T> const &cubes, float pad_dist=0.0) { // T must derive from cube_t
	cube_t tc(bcube);
	tc.expand_by_xy(pad_dist);
	for (cube_t const &c : cubes) {if (c.intersects_xy(tc)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_xy_no_adj(cube_t const &bcube, vector<T> const &cubes) { // T must derive from cube_t
	for (cube_t const &c : cubes) {if (c.intersects_xy_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> static bool check_vect_cube_contains_pt_xy(vector<T> const &cubes, point const &pos) {
	for (cube_t const &c : cubes) {if (c.contains_pt_xy(pos)) return 1;}
	return 0;
}
inline bool point_in_cubes_xy_exp(vect_cube_t const &cubes, point const &p, float dist) {
	for (cube_t const &c : cubes) {if (c.contains_pt_xy_exp(p, dist)) return 1;}
	return 0;
}
template<typename T> static bool check_vect_cube_contains_pt(vector<T> const &cubes, point const &pos) {
	for (cube_t const &c : cubes) {if (c.contains_pt(pos)) return 1;}
	return 0;
}
template<typename T> cube_t get_bcubes_union(vector<T> const &cubes) {
	if (cubes.size() == 1) {return cubes.front();}
	cube_t bcube;
	for (cube_t const &c : cubes) {bcube.assign_or_union_with_cube(c);}
	return bcube;
}
inline void swap_cube_dims(cube_t &c, unsigned d1, unsigned d2) {
	for (unsigned d = 0; d < 2; ++d) {swap(c.d[d1][d], c.d[d2][d]);}
}

struct cube_by_sz { // sort cube by size in dim descending
	bool dim;
	cube_by_sz(bool dim_) : dim(dim_) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return (b.get_sz_dim(dim) < a.get_sz_dim(dim));}
};

template<typename T> void add_to_and_clear(T &src, T &dest) {
	vector_add_to(src, dest);
	src.clear();
}
template<typename T> void add_inverted_triangles(T &verts, vector<unsigned> &indices, unsigned verts_start, unsigned ixs_start);
template<typename T> void reserve_extra(vector<T> &v, unsigned num) {v.reserve(v.size() + num);}

colorRGBA const DARK_BRASS_C(0.4, 0.35, 0.15, 1.0);

inline point get_camera_building_space() {return (get_camera_pos() - get_tiled_terrain_model_xlate());}
inline void set_cube_zvals(cube_t &c, float z1, float z2) {c.z1() = z1; c.z2() = z2;}
inline void copy_zvals(cube_t &to, cube_t const &from) {to.z1() = from.z1(); to.z2() = from.z2();}
inline float get_tc_leg_width(cube_t const &c, float width) {return min(0.5f*width*(c.dx() + c.dy()), 0.2f*c.dz());} // make legs square
inline unsigned get_rgeom_sphere_ndiv(bool low_detail) {return (low_detail ? N_SPHERE_DIV/2 : N_SPHERE_DIV);}
inline point cube_bot_center(cube_t const &c) {return point(c.xc(), c.yc(), c.z1());}
inline point cube_top_center(cube_t const &c) {return point(c.xc(), c.yc(), c.z2());}
inline bool is_known_metal_color(colorRGBA const &c) {return (c == COPPER_C || c == BRASS_C || c == DARK_BRASS_C || c == BRONZE_C || c == GOLD);}
inline colorRGBA get_specular_color(colorRGBA const &c) {return (is_known_metal_color(c) ? c : WHITE);}

void get_city_plot_zones(vect_city_zone_t &zones);
void get_city_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_city_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t const &state);
bool city_single_cube_visible_check(point const &pos, cube_t const &c);
void add_city_building_signs(cube_t const &region_bcube, vector<sign_t     > &signs);
void add_city_building_flags(cube_t const &region_bcube, vector<city_flag_t> &flags);
cube_t get_building_lights_bcube();
cube_t get_grid_bcube_for_building(building_t const &b);
unsigned get_street_dir(cube_t const &inner, cube_t const &outer);
float get_closet_wall_thickness(room_object_t const &c);
void get_closet_cubes(room_object_t const &c, cube_t cubes[5], bool for_collision=0);
void get_bed_cubes   (room_object_t const &c, cube_t cubes[6]);
void get_table_cubes (room_object_t const &c, cube_t cubes[5]);
void get_cubes_for_plastic_table(room_object_t const &c, float top_dz, cube_t cubes[3]);
void get_conf_table_cubes(room_object_t const &c, cube_t cubes[2]);
unsigned get_table_like_object_cubes(room_object_t const &c, cube_t cubes[7]);
void get_chair_cubes (room_object_t const &c, cube_t cubes[3]);
void get_tc_leg_cubes(cube_t const &c, room_object_t const &obj, float width, bool recessed, cube_t cubes[4]);
void get_reception_desk_cubes(room_object_t const &c, cube_t cubes[3]);
void get_bookcase_cubes(room_object_t const &c, cube_t &top, cube_t &middle, cube_t &back, cube_t lr[2], bool no_shelves=0, float sides_scale=1.0);
float get_drawer_cubes(room_object_t const &c, vect_cube_t &drawers, bool front_only, bool inside_only);
unsigned get_bench_cubes(room_object_t const &c, cube_t cubes[4]);
void get_diving_board_cubes(room_object_t const &c, cube_t cubes[2]);
unsigned get_shelves_for_object(room_object_t const &c, cube_t shelves[MAX_SHELVES]);
unsigned get_shelf_rack_cubes(room_object_t const &c, cube_t &back, cube_t &top, cube_t sides[2], cube_t shelves[5]);
cube_t get_shower_tub_wall(room_object_t const &c);
cube_t get_open_closet_door(room_object_t const &obj);
cube_t get_pool_table_top_surface(room_object_t const &c);
cube_t get_pan_bcube_inc_handle(room_object_t const &c);
void get_cabinet_or_counter_doors(room_object_t const &c, vect_cube_t &doors, vect_cube_t &drawers);
bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher);
room_object_t split_cabinet_at_dishwasher(room_object_t &cabinet, cube_t const &dishwasher);
room_object_t get_dresser_middle(room_object_t const &c);
room_object_t get_desk_drawers_part(room_object_t const &c);
room_object_t get_desk_top_back(room_object_t const &c);
cube_t get_attic_access_door_cube(room_object_t const &c, bool inc_ladder=0);
cube_t get_ladder_bcube_from_open_attic_door(room_object_t const &c, cube_t const &door);
cube_t get_elevator_car_panel(room_object_t const &c, float fc_thick_scale);
cube_t get_true_room_obj_bcube(room_object_t const &c);
cube_t get_sink_cube(room_object_t const &c);
cube_t get_mwave_panel_bcube(room_object_t const &c);
tquad_t get_ramp_tquad(room_object_t const &c);
colorRGBA gen_vase_color(rand_gen_t &rgen);
colorRGBA choose_pipe_color(rand_gen_t &rgen);
void gen_crate_sz(vector3d &sz, rand_gen_t &rgen, float window_vspacing);
void get_balcony_cubes(room_object_t const &c, cube_t cubes[4]);
void set_rand_pos_for_sz(cube_t &c, bool dim, float length, float width, rand_gen_t &rgen);
bool door_opens_inward(door_base_t const &door, cube_t const &room);
bool is_cube_close_to_door(cube_t const &c, float dmin, bool inc_open, cube_t const &door, unsigned check_dirs=2, unsigned open_dirs=2, bool allow_block_door=0);
void add_building_interior_lights(point const &xlate, cube_t &lights_bcube, bool sec_camera_mode);
unsigned calc_num_floors(cube_t const &c, float window_vspacing, float floor_thickness);
unsigned calc_num_floors_room(room_t const &r, float window_vspacing, float floor_thickness);
void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim);
void resize_around_center_xy(cube_t &c, float radius);
bool is_val_inside_window(cube_t const &c, bool dim, float val, float window_spacing, float window_border);
bool get_fire_ext_height_and_radius(float window_vspacing, float &height, float &radius);
point get_warning_light_src_pos(room_object_t const &c);
void offset_hanging_tv(room_object_t &obj);
template<typename T> void subtract_cube_from_cube(T const &c, cube_t const &s, vector<T> &out, bool clear_out=0);
template<typename T> void subtract_cube_from_cube_inplace(cube_t const &s, vector<T> &cubes, unsigned &ix, unsigned &iter_end);
template<typename T> void subtract_cubes_from_cube(cube_t const &c, vector<T> const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode=0);
template<typename T> bool subtract_cube_from_cubes(cube_t const &s, vector<T> &cubes, vect_cube_t *holes=nullptr, bool clip_in_z=0, bool include_adj=0, bool no_z_test=0);
template<typename T> bool line_int_cubes(point const &p1, point const &p2, vector<T> const &cubes, cube_t const &line_bcube);
void expand_to_nonzero_area(cube_t &c, float exp_amt, bool dim);
bool do_sphere_coll_polygon_sides(point &pos, cube_t const &part, float radius, bool interior_coll, vector<point> const &points, vector3d *cnorm);
int get_rect_panel_tid();
int get_bath_wind_tid ();
int get_int_door_tid  ();
int get_off_door_tid  ();
int get_bldg_door_tid ();
int get_concrete_tid  ();
int get_plywood_tid   ();
int get_insulation_tid();
int get_met_plate_tid ();
int get_mplate_nm_tid ();
int get_normal_map_for_bldg_tid(int tid);
bool has_office_chair_model();
unsigned register_sign_text(std::string const &text);
void setup_building_draw_shader(shader_t &s, float min_alpha, bool enable_indir, bool force_tsl, int use_texgen, float water_damage=0.0, float crack_damage=0.0);
void rotate_verts(vector<rgeom_mat_t::vertex_t> &verts, building_t const &building);
void add_tquad_to_verts(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
	colorRGBA const &color, vect_vnctcc_t &verts, bool invert_tc_x=0, bool exclude_frame=0, bool no_tc=0, bool no_rotate=0, bool swap_tc_xy=0);
void get_road_segs_in_region(cube_t const &region, vect_cube_t &out);
bool check_buildings_cube_coll(cube_t const &c, bool xy_only=0, bool inc_basement=1, building_t const *exclude1=nullptr, building_t const *exclude2=nullptr);
bool have_buildings_ext_paint();
void draw_buildings_ext_paint(shader_t &s);
void subtract_cube_xy(cube_t const &c, cube_t const &r, cube_t *out);
void accumulate_shared_xy_area(cube_t const &c, cube_t const &sc, float &area);
bool have_secondary_buildings();
bool get_building_door_pos_closest_to(unsigned building_id, point const &target_pos, point &door_pos, bool inc_garage_door=0);
cube_t register_deck_and_get_part_bounds(unsigned building_id, cube_t const &deck);
bool register_achievement(std::string const &str);
bool enable_building_indir_lighting_no_cib();
bool enable_building_indir_lighting();
bool in_building_gameplay_mode();
bool player_in_windowless_building();
bool player_cant_see_outside_building();
bool player_take_damage(float damage_scale, bool scream=0, int poison_type=0, uint8_t *has_key=nullptr);
void apply_building_fall_damage(float delta_z);
float get_bldg_player_height();
float get_player_eye_height();
cube_t get_stairs_bcube_expanded(stairwell_t const &s, float ends_clearance, float sides_clearance, float doorway_width);
float get_door_open_dist();
float get_lamp_width_scale();
unsigned get_face_mask(unsigned dim, bool dir);
unsigned get_skip_mask_for_xy(bool dim);
// functions in building_room_obj_expand.cc
point gen_xy_pos_in_area(cube_t const &S, vector3d const &sz, rand_gen_t &rgen, float zval=0.0);
point gen_xy_pos_in_area(cube_t const &S, float radius, rand_gen_t &rgen, float zval=0.0);
void gen_xy_pos_in_cube(point &pos, cube_t const &c, rand_gen_t &rgen);
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen, bool place_at_z1=0);
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen, bool place_at_z1=0);
// functions in building_interact.cc and building_gameplay.cc
void gen_sound_thread_safe(unsigned id, point const &pos, float gain=1.0, float pitch=1.0, float gain_scale=1.0, bool skip_if_already_playing=0);
// from building_room_geom.cpp
template<typename T> void add_sign_text_verts(std::string const &text, cube_t const &sign, bool dim, bool dir, colorRGBA const &color,
	vector<T> &verts_out, float first_char_clip_val=0.0, float last_char_clip_val=0.0, bool include_space_chars=0, bool invert_z=0);

inline void gen_sound_thread_safe_at_player(unsigned id, float gain=1.0, float pitch=1.0, bool skip_if_already_playing=0) {
	gen_sound_thread_safe(id, get_camera_pos(), gain, pitch, 1.0, skip_if_already_playing);
}
void register_building_sound(point const &pos, float volume);
void register_building_sound_at_player(float volume);
bldg_obj_type_t const &get_room_obj_type(room_object_t const &obj);
void register_building_water_splash(point const &pos, float size=1.0, bool alert_zombies=1);
bool add_water_splash(point const &pos, float radius, float size);
// functions in city_gen.cc
void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, int use_smap, int use_bmap,
	float min_alpha=0.0, bool force_tsl=0, float pcf_scale=1.0, int use_texgen=0, bool indir_lighting=0, bool is_outside=1);
void enable_animations_for_shader(shader_t &s);
void setup_city_lights(vector3d const &xlate);
void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only);
void gen_and_draw_people_in_building(ped_draw_vars_t const &pdv);
vector3d get_nom_car_size();
bool car_can_fit(cube_t const &c);
void create_mirror_reflection_if_needed(building_t const *vis_conn_bldg, vector3d const &xlate);
void draw_city_roads(int trans_op_mask, vector3d const &xlate);
void get_closest_dim_dir_xy(cube_t const &inner, cube_t const &outer, bool &dim, bool &dir);
bool check_city_tline_cube_intersect_xy(cube_t const &c);
bool connect_to_nearest_transmission_line(point const &pos, float dmax, point &conn_pt);
inline uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}
std::string gen_random_name(rand_gen_t &rgen, unsigned min_len=0, bool for_universe=0); // from Universe_name.cpp
void do_xy_rotate(float rot_sin, float rot_cos, point const &center, point &pos);
void do_xy_rotate_normal(float rot_sin, float rot_cos, point &pos);
void register_debug_event(point const &pos, std::string const &str="");

struct cmp_by_tile { // not the most efficient solution, but no memory overhead
	bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
};

