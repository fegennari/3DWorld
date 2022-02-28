// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h" // for vbo_wrap_t
#include "transform_obj.h" // for xform_matrix
#include "draw_utils.h" // for quad_batch_draw
#include "file_utils.h" // for kw_to_val_map_t

bool const EXACT_MULT_FLOOR_HEIGHT = 1;
bool const ENABLE_MIRROR_REFLECTIONS = 1;
bool const SPLIT_DOOR_PER_FLOOR    = 1; // allows mixed open/closed doors per-floor, and better texture scaling, but slower, and uses more memory
unsigned const MAX_CYLIN_SIDES     = 36;
unsigned const MAX_DRAW_BLOCKS     = 8; // for building interiors only; currently have floor, ceiling, walls, and doors
unsigned const NUM_STAIRS_PER_FLOOR= 12;
unsigned const NUM_STAIRS_PER_FLOOR_U = 16;
float const FLOOR_THICK_VAL_HOUSE  = 0.10; // 10% of floor spacing
float const FLOOR_THICK_VAL_OFFICE = 0.11; // thicker for office buildings
float const FLOOR_THICK_VAL_WINDOWLESS = 0.12; // even thicker for windowless office buildings
float const WALL_THICK_VAL         = 0.05; // 5% of floor spacing
float const DOOR_THICK_TO_WIDTH    = 0.04; // ratio of door thickness to width for doors opening to the side
float const DEF_CITY_MIN_ALPHA     = 0.01;

unsigned const NUM_BOTTLE_TYPES = 5;
unsigned const NUM_BOOK_COLORS  = 16;
unsigned const NUM_PAPER_COLORS = 6;
unsigned const NUM_SPCAN_COLORS = 10;
unsigned const NUM_LAMP_COLORS  = 6;
unsigned const NUM_TCAN_COLORS  = 6;
unsigned const NUM_TAPE_COLORS  = 7;
unsigned const NUM_SHIRT_COLORS  = 14;
colorRGBA const book_colors [NUM_BOOK_COLORS ] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, DK_BLUE, BLUE, LT_BLUE, DK_RED, RED, ORANGE, YELLOW, DK_GREEN, LT_BROWN, BROWN, DK_BROWN};
colorRGBA const spcan_colors[NUM_SPCAN_COLORS] = {WHITE, RED, GREEN, BLUE, YELLOW, PINK, ORANGE, PURPLE, BROWN, BLACK};
colorRGBA const lamp_colors[NUM_LAMP_COLORS]   = {WHITE, GRAY_BLACK, BROWN, LT_BROWN, DK_BROWN, OLIVE};
colorRGBA const cream(0.9, 0.9, 0.8), vlt_yellow(1.0, 1.0, 0.5);
colorRGBA const paper_colors[NUM_PAPER_COLORS] = {WHITE, WHITE, WHITE, cream, cream, vlt_yellow};
colorRGBA const pen_colors   [4] = {WHITE, BLACK, colorRGBA(0.2, 0.4, 1.0), RED};
colorRGBA const pencil_colors[2] = {colorRGBA(1.0, 0.75, 0.25), colorRGBA(1.0, 0.5, 0.1)};
colorRGBA const marker_colors[8] = {BLACK, RED, BLACK, BLUE, BLACK, GREEN, RED, PURPLE};
colorRGBA const tcan_colors [NUM_TCAN_COLORS ] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
colorRGBA const tape_colors [NUM_TAPE_COLORS ] = {GRAY, GRAY, GRAY, GRAY, BKGRAY, BKGRAY, colorRGBA(0.2, 0.2, 1.0)}; // gray duct tape is the most common
colorRGBA const shirt_colors[NUM_SHIRT_COLORS] = {WHITE, WHITE, WHITE, BKGRAY, BKGRAY, GRAY, GRAY, RED, BLUE, DK_BLUE, DK_GREEN, DK_BROWN, BROWN, ORANGE};
colorRGBA const LAMP_COLOR(1.0, 0.8, 0.6); // soft white
colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color; typical value to use

inline colorRGBA gen_box_color(rand_gen_t &rgen) {return colorRGBA(rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0), rgen.rand_uniform(0.9, 1.0));} // add minor color variation

class light_source;
class lmap_manager_t;
class building_nav_graph_t;
struct pedestrian_t;
struct building_t;
class building_creator_t;
class light_ix_assign_t;
struct elevator_t;
class brg_batch_draw_t;
typedef vector<vert_norm_comp_tc_color> vect_vnctcc_t;

struct bottle_params_t {
	std::string name, texture_fn;
	colorRGBA color;
	float value, label_tscale;
	bottle_params_t(std::string const &n, std::string const &fn, colorRGBA const &c, float v, float ts) : name(n), texture_fn(fn), color(c), value(v), label_tscale(ts) {}
};

// Note: we could add colorRGBA(0.8, 0.9, 1.0, 0.4) for water bottles, but transparent objects require removing interior faces such as half of the sphere
bottle_params_t const bottle_params[NUM_BOTTLE_TYPES] = {
	bottle_params_t("bottle of water",  "interiors/arrowhead_logo.jpg", colorRGBA(0.4, 0.7, 1.0 ), 1.0, 1.0),
	bottle_params_t("bottle of Coke",   "interiors/coke_label.jpg",     colorRGBA(0.2, 0.1, 0.05), 1.0, 1.0),
	bottle_params_t("bottle of beer",   "interiors/heineken_label.jpg", colorRGBA(0.1, 0.4, 0.1 ), 3.0, 2.0),
	bottle_params_t("bottle of wine",   "interiors/wine_label.jpg",     BLACK,                    10.0, 2.0),
	bottle_params_t("bottle of poison", "yuck.png",                     BLACK,                     5.0, 2.0),
};

struct rat_t {
	point pos, last_pos, dest, fear_pos;
	vector3d dir;
	float radius, height, hwidth, speed, fear, anim_time, wake_time, dist_since_sleep;
	unsigned rat_id;
	bool is_hiding, near_player, attacking;

	// this first destructor is for the lower_bound() call in vect_rat_t::get_first_rat_with_x2_gt()
	rat_t(float xval) : pos(xval, 0.0, 0.0), radius(0), height(0), hwidth(0), speed(0), fear(0),
		anim_time(0), wake_time(0), dist_since_sleep(0), rat_id(0), is_hiding(0), near_player(0), attacking(0) {}
	rat_t(point const &pos_, float radius_, vector3d const &dir_);
	bool operator<(rat_t const &r) const {return (pos.x < r.pos.x);} // compare only xvals
	bool is_moving   () const {return (speed > 0.0);}
	bool is_sleeping () const {return (wake_time > 0.0);}
	float get_hlength() const {return radius;} // this is the bounding radius, so it represents the longest dim (half length)
	point get_center () const {return point(pos.x, pos.y, (pos.z + 0.5f*height));}
	cube_t get_bcube () const; // used for collision detection and VFC; bounding cube across rotations
	cube_t get_bcube_with_dir() const; // used for model drawing; must be correct aspect ratio
	void sleep_for(float time_secs_min, float time_secs_max);
	void move(float timestep);
};

struct vect_rat_t : public vector<rat_t> {
	bool placed;
	float max_radius, max_xmove;
	vect_rat_t() : placed(0), max_radius(0.0), max_xmove(0.0) {}
	void add(rat_t const &rat);
	const_iterator get_first_rat_with_xv_gt(float x) const {return std::lower_bound(begin(), end(), rat_t(x));}
};

struct building_occlusion_state_t {
	int exclude_bix;
	bool skip_cont_camera;
	point pos;
	vector3d xlate;
	vector<cube_with_ix_t> building_ids;
	vector<point> temp_points;
	building_occlusion_state_t() : exclude_bix(-1), skip_cont_camera(0) {}

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
	bool is_occluded(cube_t const &c); // Note: non-const - state temp_points is modified
};

class occlusion_checker_noncity_t {
	building_occlusion_state_t state;
	building_creator_t const &bc;
public:
	occlusion_checker_noncity_t(building_creator_t const &bc_) : bc(bc_) {}
	void set_exclude_bix(int exclude_bix) {state.exclude_bix = exclude_bix;}
	void set_camera(pos_dir_up const &pdu);
	bool is_occluded(cube_t const &c); // Note: non-const - state temp_points is modified
};

struct city_zone_t : public cube_t {
	float zval;
	bool is_park, is_residential;
	uint8_t street_dir; // encoded as 2*dim + dir + 1; 0 is unassigned
	unsigned nbuildings, capacity; // in number of buildings; 0 is unlimited
	unsigned max_floors; // 0=unlimited
	int parent_plot_ix; // if this is a sub-plot; -1 otherwise

	city_zone_t() : zval(0.0), is_park(0), is_residential(0), street_dir(0), nbuildings(0), capacity(0), max_floors(0), parent_plot_ix(-1) {}
	city_zone_t(cube_t const &c, float zval_=0.0, bool p=0, bool r=0, unsigned sdir=0, unsigned cap=0, int ppix=-1, unsigned mf=0) :
		cube_t(c), zval(zval_), is_park(p), is_residential(r), street_dir(sdir), nbuildings(0), capacity(cap), max_floors(mf), parent_plot_ix(ppix) {}
	bool is_full() const {return (capacity > 0 && nbuildings >= capacity);}
};

typedef vector<city_zone_t> vect_city_zone_t;

struct tid_nm_pair_dstate_t {
	shader_t &s;
	int bmm_loc;
	float bump_map_mag;
	tid_nm_pair_dstate_t(shader_t &s_) : s(s_), bmm_loc(-1), bump_map_mag(1.0) {}
	void set_for_shader(float new_bump_map_mag);
	~tid_nm_pair_dstate_t();
};

struct tid_nm_pair_t { // size=28

	int tid, nm_tid; // Note: assumes each tid has only one nm_tid
	float tscale_x, tscale_y, txoff, tyoff, emissive;
	unsigned char spec_mag, shininess; // Note: spec_mag is divided by 255.0
	bool shadowed; // Note: doesn't directly affect rendering, only used for uniquing/operator==()
	bool transparent; // used to draw batched alpha blended materials last

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale_x(1.0), tscale_y(1.0), txoff(0.0), tyoff(0.0), emissive(0.0), spec_mag(0), shininess(0), shadowed(0), transparent(0) {}
	tid_nm_pair_t(int tid_, float txy=1.0, bool shadowed_=0, bool transparent_=0) : tid(tid_), nm_tid(FLAT_NMAP_TEX), tscale_x(txy), tscale_y(txy),
		txoff(0.0), tyoff(0.0), emissive(0.0), spec_mag(0), shininess(0), shadowed(shadowed_), transparent(transparent_) {} // non-normal mapped 1:1 texture AR
	tid_nm_pair_t(int tid_, int nm_tid_, float tx, float ty, float xo=0.0, float yo=0.0, bool shadowed_=0, bool transparent_=0) :
		tid(tid_), nm_tid(nm_tid_), tscale_x(tx), tscale_y(ty), txoff(xo), tyoff(yo), emissive(0.0), spec_mag(0), shininess(0), shadowed(shadowed_), transparent(transparent_) {}
	void set_specular(float mag, float shine) {spec_mag = (unsigned char)(CLIP_TO_01(mag)*255.0f); shininess = (unsigned char)max(1, min(255, round_fp(shine)));}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}

	bool is_compat_ignore_shadowed(tid_nm_pair_t const &t) const {
		return (tid == t.tid && nm_tid == t.nm_tid && emissive == t.emissive && spec_mag == t.spec_mag && shininess == t.shininess && transparent == t.transparent);
	}
	bool is_compatible(tid_nm_pair_t const &t) const {return (is_compat_ignore_shadowed(t) && shadowed == t.shadowed);}
	bool operator==(tid_nm_pair_t const &t) const {return (is_compatible(t) && tscale_x == t.tscale_x && tscale_y == t.tscale_y && txoff == t.txoff && tyoff == t.tyoff);}
	bool operator!=(tid_nm_pair_t const &t) const {return !operator==(t);}
	int get_nm_tid() const {return ((nm_tid < 0) ? FLAT_NMAP_TEX : nm_tid);}
	colorRGBA get_avg_color() const {return texture_color(tid);}
	tid_nm_pair_t get_scaled_version(float scale) const;
	bool bind_reflection_shader() const;
	void set_gl(tid_nm_pair_dstate_t &state) const;
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

	float grayscale_rand;
	colorRGBA cmin, cmax; // alpha is unused?

	color_range_t() : grayscale_rand(0.0), cmin(WHITE), cmax(WHITE) {}
	void gen_color(colorRGBA &color, rand_gen_t &rgen) const;
};

struct building_mat_t : public building_tex_params_t {

	bool no_city, add_windows, add_wind_lights;
	unsigned min_levels, max_levels, min_sides, max_sides;
	float place_radius, max_delta_z, max_rot_angle, min_level_height, min_alt, max_alt, house_prob, house_scale_min, house_scale_max;
	float split_prob, cube_prob, round_prob, asf_prob, min_fsa, max_fsa, min_asf, max_asf, wind_xscale, wind_yscale, wind_xoff, wind_yoff;
	float floor_spacing, floorplan_wind_xscale; // these are derived values
	cube_t pos_range, prev_pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color; // exterior
	colorRGBA window_color, wall_color, ceil_color, floor_color, house_ceil_color, house_floor_color;

	building_mat_t() : no_city(0), add_windows(0), add_wind_lights(0), min_levels(1), max_levels(1), min_sides(4), max_sides(4), place_radius(0.0),
		max_delta_z(0.0), max_rot_angle(0.0), min_level_height(0.0), min_alt(-1000), max_alt(1000), house_prob(0.0), house_scale_min(1.0), house_scale_max(1.0),
		split_prob(0.0), cube_prob(1.0), round_prob(0.0), asf_prob(0.0), min_fsa(0.0), max_fsa(0.0), min_asf(0.0), max_asf(0.0), wind_xscale(1.0),
		wind_yscale(1.0), wind_xoff(0.0), wind_yoff(0.0), floor_spacing(0.0), floorplan_wind_xscale(0.0), pos_range(-100,100,-100,100,0,0), prev_pos_range(all_zeros),
		sz_range(1,1,1,1,1,1), window_color(GRAY), wall_color(WHITE), ceil_color(WHITE), floor_color(LT_GRAY), house_ceil_color(WHITE), house_floor_color(WHITE) {}
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

	bool flatten_mesh, has_normal_map, tex_mirror, tex_inv_y, tt_only, infinite_buildings, dome_roof, onion_roof, enable_people_ai;
	bool gen_building_interiors, add_city_interiors, enable_rotated_room_geom, add_secondary_buildings, add_office_basements;
	unsigned num_place, num_tries, cur_prob, max_shadow_maps, buildings_rand_seed;
	float ao_factor, sec_extra_spacing, player_coll_radius_scale, interior_view_dist_scale;
	float window_width, window_height, window_xspace, window_yspace; // windows
	float wall_split_thresh, max_fp_wind_xscale, max_fp_wind_yscale, open_door_prob, locked_door_prob, basement_prob, ball_prob; // interiors
	// building AI params
	bool ai_target_player, ai_follow_player;
	unsigned ai_opens_doors; // 0=don't open doors, 1=only open if player closed door after path selection; 2=always open doors
	unsigned ai_player_vis_test; // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
	// building animal params
	unsigned num_rats_min, num_rats_max, min_attack_rats;
	float rat_speed, rat_size_min, rat_size_max;
	// gameplay state
	float player_weight_limit;
	// materials
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix, mat_gen_ix_city, mat_gen_ix_nocity, mat_gen_ix_res; // {any, city_only, non_city, residential}
	vector<unsigned> rug_tids, picture_tids, desktop_tids, sheet_tids, paper_tids;
	// use for option reading
	int read_error;
	kw_to_val_map_t<bool     >  kwmb;
	kw_to_val_map_t<unsigned >  kwmu;
	kw_to_val_map_t<float    >  kwmf;
	kw_to_val_map_t<colorRGBA>  kwmc;
	kw_to_val_map_float_check_t kwmr;

	building_params_t(unsigned num=0) : flatten_mesh(0), has_normal_map(0), tex_mirror(0), tex_inv_y(0), tt_only(0), infinite_buildings(0), dome_roof(0),
		onion_roof(0), enable_people_ai(0), gen_building_interiors(1), add_city_interiors(0), enable_rotated_room_geom(0), add_secondary_buildings(0),
		add_office_basements(0), num_place(num), num_tries(10), cur_prob(1), max_shadow_maps(32), buildings_rand_seed(0), ao_factor(0.0), sec_extra_spacing(0.0),
		player_coll_radius_scale(1.0), interior_view_dist_scale(1.0), window_width(0.0), window_height(0.0), window_xspace(0.0), window_yspace(0.0),
		wall_split_thresh(4.0), max_fp_wind_xscale(0.0), max_fp_wind_yscale(0.0), open_door_prob(1.0), locked_door_prob(0.0), basement_prob(0.5), ball_prob(0.3),
		ai_target_player(1), ai_follow_player(0), ai_opens_doors(1), ai_player_vis_test(0), num_rats_min(0), num_rats_max(0), min_attack_rats(0), rat_speed(0.0),
		rat_size_min(0.5), rat_size_max(1.0), player_weight_limit(100.0), range_translate(zero_vector), read_error(0),
		kwmb(read_error, "buildings"), kwmu(read_error, "buildings"), kwmf(read_error, "buildings"), kwmc(read_error, "buildings"), kwmr(read_error, "buildings") {init_kw_maps();}
	bool parse_buildings_option(FILE *fp);
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
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
	vector<unsigned> const &get_mat_list(bool city_only, bool non_city_only, bool residential) const {
		return (residential ? mat_gen_ix_res : (city_only ? mat_gen_ix_city : (non_city_only ? mat_gen_ix_nocity : mat_gen_ix)));
	}
	unsigned choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only, bool residential) const;
	float get_max_house_size() const;
	void set_pos_range(cube_t const &pos_range);
	void restore_prev_pos_range();
private:
	void init_kw_maps();
	int read_building_texture(FILE *fp, std::string const &str, int &error, bool check_filename=0);
	void read_texture_and_add_if_valid(FILE *fp, std::string const &str, int &error, vector<unsigned> &tids);
};

class building_draw_t;

struct building_geom_t { // describes the physical shape of a building
	unsigned num_sides;
	uint8_t door_sides[4]; // bit mask for 4 door sides, one per base part
	bool half_offset, is_pointed;
	float rot_sin, rot_cos, flat_side_amt, alt_step_factor, start_angle; // rotation in XY plane, around Z (up) axis
	//float roof_recess;

	building_geom_t(unsigned ns=4, float rs=0.0, float rc=1.0, bool pointed=0) : num_sides(ns), half_offset(0), is_pointed(pointed),
		rot_sin(rs), rot_cos(rc), flat_side_amt(0.0), alt_step_factor(0.0), start_angle(0.0)
	{
		door_sides[0] = door_sides[1] = door_sides[2] = door_sides[3] = 0;
	}
	bool is_rotated() const {return (rot_sin != 0.0);}
	bool is_cube()    const {return (num_sides == 4);}
	bool is_simple_cube()    const {return (is_cube() && !half_offset && flat_side_amt == 0.0 && alt_step_factor == 0.0);}
	bool use_cylinder_coll() const {return (num_sides > 8 && flat_side_amt == 0.0);} // use cylinder collision if not a cube, triangle, octagon, etc. (approximate)
	void do_xy_rotate    (point const &center, point &pos) const;
	void do_xy_rotate_inv(point const &center, point &pos) const;
	void do_xy_rotate_normal(point &n) const;
	void do_xy_rotate_normal_inv(point &n) const;
};

struct tquad_with_ix_t : public tquad_t {
	// roof, roof access cover, wall, chimney cap, house door, building front door, building back door, garage door, interior door back face, office door, roof door
	enum {TYPE_ROOF=0, TYPE_ROOF_ACC, TYPE_WALL, TYPE_HDOOR, TYPE_BDOOR, TYPE_BDOOR2, TYPE_GDOOR,
		TYPE_IDOOR, TYPE_IDOOR2, TYPE_ODOOR, TYPE_ODOOR2, TYPE_RDOOR, TYPE_HELIPAD, TYPE_SOLAR, TYPE_TRIM};
	bool is_building_door() const {return (type == TYPE_BDOOR || type == TYPE_BDOOR2);} // for office buildings
	bool is_exterior_door() const {return (type == TYPE_HDOOR || type == TYPE_GDOOR  || type == TYPE_RDOOR || is_building_door());}
	bool is_interior_door() const {return (type == TYPE_IDOOR || type == TYPE_IDOOR2 || type == TYPE_ODOOR || type == TYPE_ODOOR2);}

	unsigned type;
	tquad_with_ix_t(unsigned npts_=0, unsigned type_=TYPE_ROOF) : tquad_t(npts_), type(type_) {}
	tquad_with_ix_t(tquad_t const &t, unsigned type_) : tquad_t(t), type(type_) {}
};

struct vertex_range_t {
	int draw_ix; // -1 is unset
	unsigned start, end;
	vertex_range_t() : draw_ix(-1), start(0), end(0) {}
	vertex_range_t(unsigned s, unsigned e, int ix=-1) : draw_ix(ix), start(s), end(e) {}
};

struct draw_range_t {
	// intended for building interiors, which don't have many materials; may need to increase MAX_DRAW_BLOCKS later
	vertex_range_t vr[MAX_DRAW_BLOCKS]; // quad verts only for now
};

enum {
	TYPE_NONE=0, TYPE_TABLE, TYPE_CHAIR, TYPE_STAIR, TYPE_STAIR_WALL, TYPE_ELEVATOR, TYPE_LIGHT, TYPE_RUG, TYPE_PICTURE, TYPE_WBOARD,
	TYPE_BOOK, TYPE_BCASE, TYPE_TCAN, TYPE_DESK, TYPE_BED, TYPE_WINDOW, TYPE_BLOCKER, TYPE_COLLIDER, TYPE_CUBICLE, TYPE_STALL,
	TYPE_SIGN, TYPE_COUNTER, TYPE_CABINET, TYPE_KSINK, TYPE_BRSINK, TYPE_PLANT, TYPE_DRESSER, TYPE_NIGHTSTAND, TYPE_FLOORING, TYPE_CLOSET,
	TYPE_WALL_TRIM, TYPE_RAILING, TYPE_CRATE, TYPE_BOX, TYPE_MIRROR, TYPE_SHELVES, TYPE_KEYBOARD, TYPE_SHOWER, TYPE_RDESK, TYPE_BOTTLE,
	TYPE_WINE_RACK, TYPE_COMPUTER, TYPE_MWAVE, TYPE_PAPER, TYPE_BLINDS, TYPE_PEN, TYPE_PENCIL, TYPE_PAINTCAN, TYPE_LG_BALL, TYPE_HANGER_ROD,
	TYPE_DRAIN, TYPE_MONEY, TYPE_PHONE, TYPE_TPROLL, TYPE_SPRAYCAN, TYPE_MARKER, TYPE_BUTTON, TYPE_CRACK, TYPE_SWITCH, TYPE_PLATE,
	TYPE_LAPTOP, TYPE_FPLACE, TYPE_LBASKET, TYPE_WHEATER, TYPE_TAPE, TYPE_OUTLET, TYPE_PG_WALL,
	/* these next ones are all 3D models - see logic in room_object_t::is_obj_model_type() */
	TYPE_TOILET, TYPE_SINK, TYPE_TUB, TYPE_FRIDGE, TYPE_STOVE, TYPE_TV, TYPE_MONITOR, TYPE_COUCH, TYPE_OFF_CHAIR, TYPE_URINAL,
	TYPE_LAMP, TYPE_WASHER, TYPE_DRYER, TYPE_KEY, TYPE_HANGER, TYPE_CLOTHES, TYPE_FESCAPE, TYPE_CUP, TYPE_TOASTER, TYPE_RAT,
	NUM_ROBJ_TYPES};
typedef uint8_t room_object;
enum {SHAPE_CUBE=0, SHAPE_CYLIN, SHAPE_SPHERE, SHAPE_STAIRS_U, SHAPE_TALL, SHAPE_SHORT, SHAPE_ANGLED};
typedef uint8_t room_obj_shape;
enum {RTYPE_NOTSET=0, RTYPE_HALL, RTYPE_STAIRS, RTYPE_OFFICE, RTYPE_BATH, RTYPE_BED, RTYPE_KITCHEN, RTYPE_LIVING, RTYPE_DINING, RTYPE_STUDY,
	  RTYPE_ENTRY, RTYPE_LIBRARY, RTYPE_STORAGE, RTYPE_GARAGE, RTYPE_SHED, RTYPE_LOBBY, RTYPE_LAUNDRY, RTYPE_CARD, RTYPE_PLAY, RTYPE_ART,
	  RTYPE_PARKING, NUM_RTYPES};
typedef uint8_t room_type;
std::string const room_names[NUM_RTYPES] =
	{"Unset", "Hallway", "Stairs", "Office", "Bathroom", "Bedroom", "Kitchen", "Living Room", "Dining Room", "Study",
	 "Entryway", "Library", "Storage Room", "Garage", "Shed", "Lobby", "Laundry Room", "Card Room", "Play Room", "Art Room",
	 "Parking Garage"};
enum {SHAPE_STRAIGHT=0, SHAPE_U, SHAPE_WALLED, SHAPE_WALLED_SIDES};
typedef uint8_t stairs_shape;
enum {ROOM_WALL_INT=0, ROOM_WALL_SEP, ROOM_WALL_EXT, ROOM_WALL_BASEMENT};
enum {/*building models*/ OBJ_MODEL_TOILET=0, OBJ_MODEL_SINK, OBJ_MODEL_TUB, OBJ_MODEL_FRIDGE, OBJ_MODEL_STOVE, OBJ_MODEL_TV, OBJ_MODEL_MONITOR, OBJ_MODEL_COUCH,
	OBJ_MODEL_OFFICE_CHAIR, OBJ_MODEL_URINAL, OBJ_MODEL_LAMP, OBJ_MODEL_WASHER, OBJ_MODEL_DRYER, OBJ_MODEL_KEY, OBJ_MODEL_HANGER, OBJ_MODEL_CLOTHES,
	OBJ_MODEL_FESCAPE, OBJ_MODEL_CUP, OBJ_MODEL_TOASTER, OBJ_MODEL_RAT,
	/*city models*/ OBJ_MODEL_FHYDRANT, OBJ_MODEL_SUBSTATION, OBJ_MODEL_UMBRELLA, NUM_OBJ_MODELS};

// object flags
unsigned const RO_FLAG_LIT     = 0x01; // light is on
unsigned const RO_FLAG_TOS     = 0x02; // at top of stairs
unsigned const RO_FLAG_RSTAIRS = 0x04; // in a room with stairs
unsigned const RO_FLAG_INVIS   = 0x08; // invisible
unsigned const RO_FLAG_NOCOLL  = 0x10; // no collision detection
unsigned const RO_FLAG_OPEN    = 0x20; // open, for elevators, closet doors, bathroom stalls, and phones
unsigned const RO_FLAG_NODYNAM = 0x40; // for light shadow maps
unsigned const RO_FLAG_INTERIOR= 0x80; // applies to containing room
// object flags, second byte
unsigned const RO_FLAG_EMISSIVE= 0x0100; // for signs, lights, and phones
unsigned const RO_FLAG_HANGING = 0x0200; // for signs, blinds, hangers, and shirts; treated as "folding" for closet doors
unsigned const RO_FLAG_ADJ_LO  = 0x0400; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_HI  = 0x0800; // for kitchen counters/closets/door trim/blinds/railings
unsigned const RO_FLAG_ADJ_BOT = 0x1000; // for door trim
unsigned const RO_FLAG_ADJ_TOP = 0x2000; // for door trim/railings
unsigned const RO_FLAG_IS_HOUSE= 0x4000; // used for mirror reflections, shelves, and tables
unsigned const RO_FLAG_RAND_ROT= 0x8000; // random rotation; used for office chairs, papers, pictures, and cups
unsigned const RO_FLAG_UNTEXTURED = 0x1000; // for shirts, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_FROM_SET   = 0x1000; // for books, aliased with RO_FLAG_ADJ_BOT
unsigned const RO_FLAG_HAS_VOL_IX = 0x2000; // for books, aliased with RO_FLAG_ADJ_TOP
// object flags, third byte, for pickup/interact state
unsigned const RO_FLAG_TAKEN1  = 0x010000; // no picture / no bed pillows
unsigned const RO_FLAG_TAKEN2  = 0x020000; // no bed sheets
unsigned const RO_FLAG_TAKEN3  = 0x040000; // no bed mattress
unsigned const RO_FLAG_TAKEN4  = 0x080000; // for future use
unsigned const RO_FLAG_EXPANDED= 0x100000; // for shelves and closets
unsigned const RO_FLAG_WAS_EXP = 0x200000; // for objects in/on shelves, closets, and drawers, cabinets, and books
unsigned const RO_FLAG_ROTATING= 0x400000; // for office chairs
unsigned const RO_FLAG_IN_CLOSET=0x800000; // for closet lights
// object flags, fourth byte
unsigned const RO_FLAG_DYNAMIC  = 0x01000000; // dynamic object (balls, elevators, etc.)
unsigned const RO_FLAG_DSTATE   = 0x02000000; // this object has dynamic state
unsigned const RO_FLAG_NO_CONS  = 0x04000000; // this object is not consumable (bottles)
unsigned const RO_FLAG_IS_ACTIVE= 0x08000000; // active, for sinks, tubs, buttons, etc.
unsigned const RO_FLAG_USED     = 0x10000000; // used by the player (spraypaint, marker, etc.)
unsigned const RO_FLAG_IN_ELEV  = 0x20000000; // for elevator lights and buttons
unsigned const RO_FLAG_BROKEN   = 0x40000000; // for TVs and monitors, maybe can use for windows
unsigned const RO_FLAG_MOVED    = 0x80000000; // for player push/pull

struct bldg_obj_type_t {
	bool player_coll=0, ai_coll=0, rat_coll=0, pickup=0, attached=0, is_model=0;
	uint8_t lg_sm=0; // 0=neither (model), 1=large item, 2=small item, 3=split into large and small
	float value=0.0, weight=0.0;
	unsigned capacity=0; // for consumable/usable objects
	std::string name;

	bldg_obj_type_t() {}
	bldg_obj_type_t(bool pc, bool ac, bool rc, bool pu, bool at, bool im, uint8_t ls, float v, float w, std::string const &n, unsigned cap=0) :
		player_coll(pc), ai_coll(ac), rat_coll(rc), pickup(pu), attached(at), is_model(im), lg_sm(ls), value(v), weight(w), capacity(cap), name(n) {}
};

struct room_object_t : public cube_t {
	bool dim, dir;
	uint8_t room_id; // for at most 256 rooms per floor
	uint16_t obj_id, drawer_flags, item_flags;
	room_object type; // 8-bit
	room_obj_shape shape; // 8-bit
	unsigned flags;
	float light_amt;
	colorRGBA color;

	room_object_t() : dim(0), dir(0), room_id(0), obj_id(0), drawer_flags(0), item_flags(0), type(TYPE_NONE), shape(SHAPE_CUBE), flags(0), light_amt(1.0) {}
	room_object_t(cube_t const &c, room_object type_, uint8_t rid, bool dim_=0, bool dir_=0, unsigned f=0, float light=1.0,
		room_obj_shape shape_=SHAPE_CUBE, colorRGBA const color_=WHITE) :
		cube_t(c), dim(dim_), dir(dir_), room_id(rid), obj_id(0), drawer_flags(0), item_flags(0), type(type_), shape(shape_), flags(f), light_amt(light), color(color_)
	{check_normalized();}
	void check_normalized() const;
	unsigned get_combined_flags() const {return (((unsigned)drawer_flags << 16) + (unsigned)item_flags);} // treat {drawer_flags, item_flags} as a single 32-bit flags
	void set_combined_flags(unsigned v) {drawer_flags = (v >> 16); item_flags = (v & 0xFFFF);}
	bool is_valid   () const {return  (type != TYPE_NONE);}
	bool is_lit     () const {return  (flags & RO_FLAG_LIT);}
	bool has_stairs () const {return  (flags & (RO_FLAG_TOS | RO_FLAG_RSTAIRS));}
	bool is_visible () const {return !(flags & RO_FLAG_INVIS);}
	bool no_coll    () const {return  (flags & RO_FLAG_NOCOLL);}
	bool is_interior() const {return  (flags & RO_FLAG_INTERIOR);}
	bool is_open    () const {return  (flags & RO_FLAG_OPEN);}
	bool is_house   () const {return  (flags & RO_FLAG_IS_HOUSE);}
	bool was_expanded() const{return  (flags & RO_FLAG_WAS_EXP);}
	bool is_dynamic () const {return  (flags & RO_FLAG_DYNAMIC);}
	bool has_dstate () const {return  (flags & RO_FLAG_DSTATE);}
	bool is_moving  () const {return (is_dynamic() && has_dstate());}
	bool was_moved  () const {return  (flags & RO_FLAG_MOVED);}
	bool is_light_type() const {return (type == TYPE_LIGHT || (type == TYPE_LAMP && !was_expanded()));} // light, or lamp not in closet
	bool is_sink_type () const {return (type == TYPE_SINK || type == TYPE_KSINK || type == TYPE_BRSINK);}
	bool is_obj_model_type() const {return (type >= TYPE_TOILET && type < NUM_ROBJ_TYPES);}
	bool is_small_closet() const {return (get_sz_dim(!dim) < 1.2*dz());}
	bool is_bottle_empty() const {return ((obj_id & 192) == 192);} // empty if both bits 6 and 7 are set
	bool desk_has_drawers()const {return bool(room_id & 3);} // 75% of the time
	bool can_use        () const;
	bool is_interactive () const {return (has_dstate() || can_use());}
	bool can_place_onto () const;
	bool is_floor_collidable() const;
	unsigned get_bottle_type() const {return ((obj_id&63) % NUM_BOTTLE_TYPES);} // first 6 bits are bottle type
	unsigned get_orient () const {return (2*dim + dir);}
	float get_radius() const;
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
	int get_model_id () const;
	void set_as_bottle(unsigned rand_id, unsigned max_type=NUM_BOTTLE_TYPES-1, bool no_empty=0);
	colorRGBA get_color() const;
	colorRGBA get_model_color() const;
	vector3d get_dir() const {vector3d v(zero_vector); v[dim] = (dir ? 1.0 : -1.0); return v;}
	void set_rand_gen_state(rand_gen_t &rgen) const {rgen.set_state(obj_id+1, room_id+1);}
};
typedef vector<room_object_t> vect_room_object_t;

struct carried_item_t : public room_object_t {
	unsigned use_count;
	carried_item_t() : use_count(0) {}
	carried_item_t(room_object_t const &o) : room_object_t(o), use_count(0) {}
	float get_remaining_capacity_ratio() const;
};

struct room_obj_dstate_t { // state used for dynamic room objects
	vector3d velocity;
	xform_matrix rot_matrix;
};

struct rgeom_storage_t {
	typedef vert_norm_comp_tc_color vertex_t;
	vector<vertex_t> quad_verts, itri_verts;
	vector<unsigned> indices; // for indexed_quad_verts
	tid_nm_pair_t tex; // used sort of like a map key

	rgeom_storage_t() {}
	rgeom_storage_t(tid_nm_pair_t const &tex_) : tex(tex_) {}
	bool empty() const {return (quad_verts.empty() && itri_verts.empty());}
	void clear();
	void swap_vectors(rgeom_storage_t &s);
	void swap(rgeom_storage_t &s);
	unsigned get_tot_vert_capacity() const {return (quad_verts.capacity() + itri_verts.capacity());}
	unsigned get_mem_usage() const {return (get_cont_mem_usage(quad_verts) + get_cont_mem_usage(itri_verts) + get_cont_mem_usage(indices));}
};

class rgeom_mat_t : public rgeom_storage_t { // simplified version of building_draw_t::draw_block_t

	indexed_vao_manager_with_shadow_t vao_mgr;
public:
	unsigned num_verts, num_ixs; // for drawing
	uint8_t dir_mask; // {-x, +x, -y, +y, -z, +z}
	bool en_shadows;

	rgeom_mat_t(tid_nm_pair_t const &tex_=tid_nm_pair_t()) : rgeom_storage_t(tex_), num_verts(0), num_ixs(0), dir_mask(0), en_shadows(0) {}
	//~rgeom_mat_t() {assert(vbo_mgr.vbo == 0); assert(vbo_mgr.ivbo == 0);} // VBOs should be freed before destruction
	void enable_shadows() {en_shadows = 1;}
	void clear();
	void add_cube_to_verts(cube_t const &c, colorRGBA const &color, point const &tex_origin=all_zeros,
		unsigned skip_faces=0, bool swap_tex_st=0, bool mirror_x=0, bool mirror_y=0, bool inverted=0);
	void add_cube_to_verts_untextured(cube_t const &c, colorRGBA const &color, unsigned skip_faces=0);
	void add_ortho_cylin_to_verts(cube_t const &c, colorRGBA const &color, int dim, bool draw_bot, bool draw_top, bool two_sided=0, bool inv_tb=0,
		float rs_bot=1.0, float rs_top=1.0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0, bool swap_txy=0);
	void add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided=0, bool inv_tb=0,
		float rs_bot=1.0, float rs_top=1.0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0, bool swap_txy=0);
	void add_cylin_to_verts(point const &bot, point const &top, float bot_radius, float top_radius, colorRGBA const &color, bool draw_bot, bool draw_top,
		bool two_sided=0, bool inv_tb=0, float side_tscale=1.0, float end_tscale=1.0, bool skip_sides=0, unsigned ndiv=N_CYL_SIDES, float side_tscale_add=0.0, bool swap_txy=0);
	void add_disk_to_verts(point const &pos, float radius, bool normal_z_neg, colorRGBA const &color);
	void add_sphere_to_verts(cube_t const &c, colorRGBA const &color, bool low_detail=0, vector3d const &skip_hemi_dir=zero_vector, xform_matrix const *const matrix=nullptr);
	void add_triangle_to_verts(point const v[3], colorRGBA const &color, bool two_sided, float tscale=1.0);
	void create_vbo(building_t const &building);
	void create_vbo_inner();
	void draw(tid_nm_pair_dstate_t &state, brg_batch_draw_t *bbd, int shadow_only, bool reflection_pass);
	void draw_inner(tid_nm_pair_dstate_t &state, int shadow_only) const;
	void upload_draw_and_clear(tid_nm_pair_dstate_t &state);
};

struct building_materials_t : public vector<rgeom_mat_t> {
	void clear();
	unsigned count_all_verts() const;
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows);
	void create_vbos(building_t const &building);
	void draw(brg_batch_draw_t *bbd, shader_t &s, int shadow_only, bool reflection_pass);
	void upload_draw_and_clear(shader_t &s);
};

class brg_batch_draw_t {
	struct mat_entry_t {
		tid_nm_pair_t tex;
		vector<rgeom_mat_t const *> mats;
		mat_entry_t() {}
		mat_entry_t(rgeom_mat_t const &m) : tex(m.tex) {mats.push_back(&m);}
	};
	vector<mat_entry_t> to_draw;
	vector<int> tid_to_first_mat_map; // -1 is unset
public:
	uint8_t camera_dir_mask;

	brg_batch_draw_t() : camera_dir_mask(0) {}
	void clear() {to_draw.clear(); tid_to_first_mat_map.clear();}
	void set_camera_dir_mask(point const &camera_bs, cube_t const &bcube);
	void add_material(rgeom_mat_t const &m);
	void draw_and_clear(shader_t &s);
};

struct obj_model_inst_t {
	unsigned obj_id;
	vector3d dir;
	obj_model_inst_t(unsigned oid, vector3d const &d) : obj_id(oid), dir(d) {}
};

struct tape_quad_batch_draw : public quad_batch_draw {
	int moving_vert_cyilin_int_tape(point &cur_pos, point const &prev_pos, float z1, float z2, float radius, float slow_amt, bool is_player) const;
	void split_tape_at(unsigned first_vert, point const &pos, float min_zval);
};

struct paint_draw_t {
	quad_batch_draw qbd[2]; // {spraypaint, markers}
	void draw_paint() const;
};
struct building_decal_manager_t {
	paint_draw_t paint_draw[2]; // {interior, exterior}
	quad_batch_draw blood_qbd, tp_qbd, pend_tape_qbd;
	tape_quad_batch_draw tape_qbd; // for tape, but not pend_tape because it hasn't been placed yet

	void commit_pend_tape_qbd();
	void draw_building_interior_decals(bool player_in_building, bool shadow_only) const;
};

struct building_room_geom_t {

	bool has_elevators, has_pictures, lights_changed, lighting_invalid, modified_by_player;
	unsigned char num_pic_tids;
	float obj_scale;
	unsigned buttons_start, stairs_start; // index of first object of {TYPE_TRIM, TYPE_BUTTON, TYPE_STAIR}
	point tex_origin;
	colorRGBA wood_color;
	// objects in rooms; expanded_objs is for things that have been expanded for player interaction; model_objs is for models in drawers; trim_objs is for wall/door/window trim
	vect_room_object_t objs, expanded_objs, model_objs, trim_objs;
	vector<room_obj_dstate_t> obj_dstate;
	vector<obj_model_inst_t> obj_model_insts;
	vector<unsigned> moved_obj_ids;
	vect_rat_t rats;
	// {large static, small static, dynamic, lights, alpha mask, transparent, door} materials
	building_materials_t mats_static, mats_small, mats_detail, mats_dynamic, mats_lights, mats_amask, mats_alpha, mats_doors;
	vect_cube_t light_bcubes;
	building_decal_manager_t decal_manager;

	building_room_geom_t(point const &tex_origin_=all_zeros) : has_elevators(0), has_pictures(0), lights_changed(0), lighting_invalid(0), modified_by_player(0),
		num_pic_tids(0), obj_scale(1.0), buttons_start(0), stairs_start(0), tex_origin(tex_origin_), wood_color(WHITE) {}
	bool empty() const {return objs.empty();}
	void clear();
	void clear_materials();
	void clear_lit_materials();
	void clear_static_vbos();
	void clear_static_small_vbos();
	void clear_and_recreate_lights() {lights_changed = 1;} // cache the state and apply the change later in case this is called from a different thread
	unsigned get_num_verts() const {return (mats_static.count_all_verts() + mats_small.count_all_verts() + mats_detail.count_all_verts() + mats_dynamic.count_all_verts() +
		mats_lights.count_all_verts() + mats_amask.count_all_verts() + mats_alpha.count_all_verts() + mats_doors.count_all_verts());}
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool transparent=0);
	rgeom_mat_t &get_untextured_material(bool inc_shadows=0, bool dynamic=0, unsigned small=0, bool transparent=0) {
		return get_material(tid_nm_pair_t(-1, 1.0, inc_shadows, transparent), inc_shadows, dynamic, small, transparent);
	}
	rgeom_mat_t &get_wood_material(float tscale=1.0, bool inc_shadows=1, bool dynamic=0, unsigned small=0);
	rgeom_mat_t &get_metal_material(bool inc_shadows=0, bool dynamic=0, unsigned small=0);
	colorRGBA apply_wood_light_color(room_object_t const &o) const;
	void add_tquad(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
		colorRGBA const &color, bool invert_tc_x, bool exclude_frame, bool no_tc);
	vect_room_object_t::const_iterator get_placed_objs_end() const {return (objs.begin() + buttons_start);} // excludes buttons, stairs, and elevators
	vect_room_object_t::const_iterator get_stairs_start   () const {return (objs.begin() + stairs_start );} // excludes stairs
	// Note: these functions are all for drawing objects / adding them to the vertex list
	void add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale, bool use_metal_mat=0, bool draw_tops=0);
	void add_table(room_object_t const &c, float tscale, float top_dz, float leg_width);
	void add_chair(room_object_t const &c, float tscale);
	void add_dresser(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void add_dresser_drawers(room_object_t const &c, float tscale);
	void add_stair(room_object_t const &c, float tscale, vector3d const &tex_origin);
	void add_interior_wall(room_object_t const &c, vector3d const &tex_origin, tid_nm_pair_t const &wall_tex, bool draw_top_bot);
	void add_elevator(room_object_t const &c, float tscale, float fc_thick_scale);
	void add_elevator_doors(elevator_t const &e, float fc_thick_scale);
	void add_light(room_object_t const &c, float tscale);
	void add_rug(room_object_t const &c);
	void add_picture(room_object_t const &c);
	void add_book_title(std::string const &title, cube_t const &title_area, rgeom_mat_t &mat, colorRGBA const &color,
		unsigned hdim, unsigned tdim, unsigned wdim, bool cdir, bool ldir, bool wdir);
	void add_book(room_object_t const &c, bool inc_lg, bool inc_sm, float tilt_angle=0.0, unsigned extra_skip_faces=0, bool no_title=0, float z_rot_angle=0.0);
	void add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale, bool no_shelves=0, float sides_scale=1.0,
		point const *const use_this_tex_origin=nullptr, vect_room_object_t *books=nullptr);
	void add_wine_rack(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale);
	void add_desk(room_object_t const &c, float tscale, bool inc_lg, bool inc_sm);
	void add_reception_desk(room_object_t const &c, float tscale);
	void add_bed(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale);
	void add_window(room_object_t const &c, float tscale);
	void add_crack(room_object_t const &c);
	void add_switch(room_object_t const &c, bool draw_detail_pass);
	void add_outlet(room_object_t const &c);
	void add_plate(room_object_t const &c);
	void add_tub_outer(room_object_t const &c);
	void add_tv_picture(room_object_t const &c);
	void add_trashcan(room_object_t const &c);
	void add_laundry_basket(room_object_t const &c);
	void add_water_heater(room_object_t const &c);
	void add_toaster_proxy(room_object_t const &c);
	void add_br_stall(room_object_t const &c);
	void add_cubicle(room_object_t const &c, float tscale);
	void add_sign(room_object_t const &c, bool inc_back, bool inc_text);
	void add_counter(room_object_t const &c, float tscale);
	void add_cabinet(room_object_t const &c, float tscale);
	void add_closet(room_object_t const &c, tid_nm_pair_t const &wall_tex, bool inc_lg, bool inc_sm);
	void add_hanger_rod(room_object_t const &c);
	void add_drain_pipe(room_object_t const &c);
	void add_key(room_object_t const &c);
	void add_money(room_object_t const &c);
	void add_phone(room_object_t const &c);
	void add_tproll(room_object_t const &c);
	void add_tape(room_object_t const &c);
	static void add_spraycan_to_material(room_object_t const &c, rgeom_mat_t &mat);
	void add_spraycan(room_object_t const &c);
	void add_button(room_object_t const &c);
	void add_crate(room_object_t const &c);
	void add_box(room_object_t const &c);
	void add_paint_can(room_object_t const &c);
	void add_shelves(room_object_t const &c, float tscale);
	void add_keyboard(room_object_t const &c);
	void add_obj_with_top_texture  (room_object_t const &c, std::string const &texture_name, colorRGBA const &sides_color, bool is_small=0);
	void add_obj_with_front_texture(room_object_t const &c, std::string const &texture_name, colorRGBA const &sides_color, bool is_small=0);
	void add_computer(room_object_t const &c);
	void add_laptop(room_object_t const &c);
	void add_mwave(room_object_t const &c);
	void add_mirror(room_object_t const &c);
	void add_shower(room_object_t const &c, float tscale);
	void add_bottle(room_object_t const &c);
	void add_paper(room_object_t const &c);
	static void add_pen_pencil_marker_to_material(room_object_t const &c_, rgeom_mat_t &mat);
	void add_pen_pencil_marker(room_object_t const &c);
	void add_flooring(room_object_t const &c, float tscale);
	void add_wall_trim(room_object_t const &c, bool for_closet=0);
	void add_blinds(room_object_t const &c);
	void add_fireplace(room_object_t const &c, float tscale);
	void add_railing(room_object_t const &c);
	void add_potted_plant(room_object_t const &c, bool inc_pot, bool inc_plant);
	void add_lg_ball(room_object_t const &c);
	static void draw_lg_ball_in_building   (room_object_t  const &c, shader_t &s);
	static void draw_interactive_player_obj(carried_item_t const &c, shader_t &s, vector3d const &xlate);
	// functions for expanding nested objects
	void expand_shelves(room_object_t const &c);
	void get_bookcase_books(room_object_t const &c, vect_room_object_t &books) {add_bookcase(c, 0, 0, 1.0, 0, 1.0, nullptr, &books);} // Note: technically const
	void expand_closet(room_object_t const &c) {add_closet_objects(c, expanded_objs);}
	void expand_cabinet(room_object_t const &c);
	void expand_wine_rack(room_object_t const &c) {add_wine_rack_bottles(c, expanded_objs);}
	void expand_object(room_object_t &c);
	static room_object_t get_item_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix);
	// other functions
	bool closet_light_is_on(cube_t const &closet) const;
	int find_nearest_pickup_object(building_t const &building, point const &at_pos, vector3d const &in_dir, float range, float &obj_dist) const;
	bool cube_intersects_moved_obj(cube_t const &c, int ignore_obj_id=-1) const;
	bool open_nearest_drawer(building_t &building, point const &at_pos, vector3d const &in_dir, float range, bool pickup_item, bool check_only);
	void remove_object(unsigned obj_id, building_t &building);
	bool player_pickup_object(building_t &building, point const &at_pos, vector3d const &in_dir);
	void update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken);
	void update_draw_state_for_obj_type_flags(bldg_obj_type_t const &type, building_t &building);
	room_object_t &get_room_object_by_index(unsigned obj_id);
	int find_avail_obj_slot() const;
	void add_expanded_object(room_object_t const &obj);
	bool add_room_object(room_object_t const &obj, building_t &building, bool set_obj_id=0, vector3d const &velocity=zero_vector);
	void update_dynamic_draw_data() {mats_dynamic.clear();}
	void create_static_vbos(building_t const &building);
	void create_small_static_vbos(building_t const &building);
	void create_detail_vbos(building_t const &building);
	void add_small_static_objs_to_verts(vect_room_object_t const &objs_to_add);
	void create_obj_model_insts(building_t const &building);
	void create_lights_vbos(building_t const &building);
	void create_dynamic_vbos(building_t const &building);
	void create_door_vbos(building_t const &building);
	void draw(brg_batch_draw_t *bbd, shader_t &s, building_t const &building, occlusion_checker_noncity_t &oc, vector3d const &xlate,
		unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building);
	unsigned allocate_dynamic_state();
	room_obj_dstate_t &get_dstate(room_object_t const &obj);
private:
	static void add_closet_objects(room_object_t const &c, vect_room_object_t &objects);
	static unsigned get_shelves_for_object(room_object_t const &c, cube_t shelves[4]);
	static void get_shelf_objects(room_object_t const &c_in, cube_t const shelves[4], unsigned num_shelves, vect_room_object_t &objects);
	static void add_wine_rack_bottles(room_object_t const &c, vect_room_object_t &objects);
	static void add_vert_roll_to_material(room_object_t const &c, rgeom_mat_t &mat, float sz_ratio=1.0, bool player_held=0);
	void add_bcase_book(room_object_t const &c, cube_t const &book, bool inc_lg, bool inc_sm, bool backwards, bool in_set,
		unsigned skip_faces, unsigned book_ix, unsigned set_start_ix, colorRGBA const &color, vect_room_object_t *books);
}; // building_room_geom_t

struct elevator_t : public cube_t {
	bool dim, dir, at_edge, was_called; // door dim/dir
	unsigned room_id, car_obj_id, light_obj_id, button_id_start, button_id_end;
	float target_zval, open_amt;

	elevator_t(cube_t const &c, unsigned rid, bool dim_, bool dir_, bool at_edge_) :
		cube_t(c), dim(dim_), dir(dir_), at_edge(at_edge_), was_called(0), room_id(rid), car_obj_id(0),
		light_obj_id(0), button_id_start(0), button_id_end(0), target_zval(0.0), open_amt(0.0)
	{assert(is_strictly_normalized());}
	float get_wall_thickness () const {return 0.02*get_sz_dim(!dim);}
	float get_frame_width    () const {return 0.20*get_sz_dim(!dim);}
	unsigned get_door_face_id() const {return (2*dim + dir);}
	unsigned get_coll_cubes(cube_t cubes[5]) const; // returns 1 or 5 cubes
	void call_elevator(float targ_z) {target_zval = targ_z; was_called = 1;}
};

unsigned const NUM_RTYPE_SLOTS = 6; // enough for houses

struct room_t : public cube_t {
	uint8_t has_stairs; // per-floor bit mask; always set to 255 for stairs that span the entire room
	bool has_elevator, has_center_stairs, no_geom, is_hallway, is_office, is_sec_bldg, interior;
	uint8_t ext_sides; // sides that have exteriors, and likely windows (bits for x1, x2, y1, y2)
	uint8_t part_id, num_lights;
	room_type rtype[NUM_RTYPE_SLOTS]; // this applies to the first floor because some rooms can have variable per-floor assignment
	uint64_t lit_by_floor;
	float light_intensity; // due to room lights, if turned on

	room_t() : has_stairs(0), has_elevator(0), has_center_stairs(0), no_geom(0), is_hallway(0), is_office(0), is_sec_bldg(0), interior(0),
		ext_sides(0), part_id(0), num_lights(0), lit_by_floor(0), light_intensity(0.0) {assign_all_to(RTYPE_NOTSET);}
	room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_);
	void assign_all_to(room_type rt);
	void assign_to(room_type rt, unsigned floor);
	room_type get_room_type(unsigned floor) const {return rtype[min(floor, NUM_RTYPE_SLOTS-1U)];}
	bool has_bathroom() const;
	bool is_lit_on_floor    (unsigned floor) const {return (lit_by_floor & (1ULL << (floor&63)));}
	bool has_stairs_on_floor(unsigned floor) const {return (has_stairs & (1U << min(floor, 7U)));} // floors >= 7 are treated as the top floor
	bool is_garage_or_shed  (unsigned floor) const {return (is_sec_bldg || get_room_type(floor) == RTYPE_GARAGE || get_room_type(floor) == RTYPE_SHED);}
	float get_light_amt() const;
	unsigned get_floor_containing_zval(float zval, float floor_spacing) const {return (is_sec_bldg ? 0 : unsigned((zval - z1())/floor_spacing));}
};

struct stairs_landing_base_t {
	bool dim, dir, roof_access, stack_conn, against_wall[2];
	stairs_shape shape;

	stairs_landing_base_t() : dim(0), dir(0), roof_access(0), stack_conn(0), shape(SHAPE_STRAIGHT) {against_wall[0] = against_wall[1] = 0;}
	stairs_landing_base_t(bool dim_, bool dir_, bool roof_access_, stairs_shape shape_, bool sc=0) :
		dim(dim_), dir(dir_), roof_access(roof_access_), stack_conn(sc), shape(shape_) {against_wall[0] = against_wall[1] = 0;}
	void set_against_wall(bool val[2]) {against_wall[0] = val[0]; against_wall[1] = val[1];}
	unsigned get_face_id   () const {return (2*dim + dir);}
	unsigned get_num_stairs() const {return ((shape == SHAPE_U) ? NUM_STAIRS_PER_FLOOR_U : NUM_STAIRS_PER_FLOOR);}
};

struct landing_t : public cube_t, public stairs_landing_base_t {
	bool for_elevator, has_railing, is_at_top;
	uint8_t floor;

	landing_t() : for_elevator(0), has_railing(0), is_at_top(0), floor(0) {}
	landing_t(cube_t const &c, bool e, uint8_t f, bool dim_, bool dir_, bool railing=0, stairs_shape shape_=SHAPE_STRAIGHT, bool roof_access_=0, bool at_top=0, bool sc=0) :
		cube_t(c), stairs_landing_base_t(dim_, dir_, roof_access_, shape_, sc), for_elevator(e), has_railing(railing), is_at_top(at_top), floor(f)
	{assert(is_strictly_normalized());}
};

struct stairwell_t : public cube_t, public stairs_landing_base_t {
	uint8_t num_floors;
	int16_t stairs_door_ix;

	stairwell_t() : num_floors(0), stairs_door_ix(-1) {}
	stairwell_t(cube_t const &c, unsigned n, bool dim_, bool dir_, stairs_shape s=SHAPE_STRAIGHT, bool r=0, bool sc=0) :
		cube_t(c), stairs_landing_base_t(dim_, dir_, r, s, sc), num_floors(n), stairs_door_ix(-1) {}
};
typedef vector<stairwell_t> vect_stairwell_t;

struct door_base_t : public cube_t {
	bool dim, open_dir, hinge_side, on_stairs;
	// is it useful to store the two rooms in the door/door_stack? this will speed up connectivity searches for navigation and room assignment,
	// but only for finding the second room connected to a door, because we still need to iterate over all doors;
	// unfortunately, it's not easy/cheap to assign these values because the room may not even be added until after the door is placed, so we have to go back and set room1/room2 later
	//uint8_t room1, room2;
	door_base_t() : dim(0), open_dir(0), hinge_side(0), on_stairs(0) {}
	door_base_t(cube_t const &c, bool dim_, bool dir, bool os=0, bool hs=0) :
		cube_t(c), dim(dim_), open_dir(dir), hinge_side(hs), on_stairs(os) {assert(is_strictly_normalized());}
	bool get_check_dirs  () const {return (dim ^ open_dir ^ hinge_side ^ 1);}
	float get_width      () const {return get_sz_dim(!dim);}
	float get_thickness  () const {return DOOR_THICK_TO_WIDTH*get_width();}
	cube_t get_true_bcube() const {cube_t bc(*this); bc.expand_in_dim(dim, 0.5*get_thickness()); return bc;}
	bool is_same_stack(door_base_t const &d) const {return (d.x1() == x1() && d.y1() == y1());}
};
struct door_stack_t : public door_base_t {
	unsigned first_door_ix; // on the lowest floor
	door_stack_t() : first_door_ix(0) {}
	door_stack_t(door_base_t const &db, unsigned fdix) : door_base_t(db), first_door_ix(fdix) {}
};
struct door_t : public door_base_t {
	bool open, locked, blocked;
	door_t() : open(0), locked(0), blocked(0) {}
	door_t(cube_t const &c, bool dim_, bool dir, bool open_=1, bool os=0, bool hs=0) : door_base_t(c, dim_, dir, os, hs), open(open_), locked(0), blocked(0) {}
	bool is_closed_and_locked() const {return (!open && locked);}
	bool is_locked_or_blocked(bool have_key) const {return (blocked || (is_closed_and_locked() && !have_key));}
};
typedef vector<door_stack_t> vect_door_stack_t;
typedef vector<door_t> vect_door_t;

enum {ROOF_OBJ_BLOCK=0, ROOF_OBJ_ANT, ROOF_OBJ_WALL, ROOF_OBJ_ECAP, ROOF_OBJ_AC, ROOF_OBJ_SCAP};
enum {ROOF_TYPE_FLAT=0, ROOF_TYPE_SLOPE, ROOF_TYPE_PEAK, ROOF_TYPE_DOME, ROOF_TYPE_ONION};

struct roof_obj_t : public cube_t {
	uint8_t type;
	roof_obj_t(uint8_t type_=ROOF_OBJ_BLOCK) : type(type_) {}
	roof_obj_t(cube_t const &c, uint8_t type_=ROOF_OBJ_BLOCK) : cube_t(c), type(type_) {assert(is_strictly_normalized());}
};
typedef vector<roof_obj_t> vect_roof_obj_t;


// may as well make this its own class, since it could get large and it won't be used for every building
struct building_interior_t {
	vect_cube_t floors, ceilings, walls[2]; // walls are split by dim
	vect_stairwell_t stairwells;
	vector<door_t> doors;
	vector<door_stack_t> door_stacks;
	vector<landing_t> landings; // for stairs and elevators
	vector<room_t> rooms;
	vector<elevator_t> elevators;
	vect_cube_t exclusion;
	std::unique_ptr<building_room_geom_t> room_geom;
	std::unique_ptr<building_nav_graph_t> nav_graph;
	draw_range_t draw_range;
	uint64_t top_ceilings_mask; // bit mask for ceilings that are on the top floor and have no floor above them
	bool door_state_updated, is_unconnected;

	building_interior_t();
	~building_interior_t();
	float get_doorway_width() const;
	bool is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin=0.0f, bool inc_open=0, bool check_open_dir=0) const;
	bool is_blocked_by_stairs_or_elevator(cube_t const &c, float dmin=0.0f, bool elevators_only=0) const;
	bool is_blocked_by_stairs_or_elevator_no_expand(cube_t const &c, float dmin=0.0f) const;
	void finalize();
	bool update_elevators(building_t const &building, point const &player_pos);
	bool check_sphere_coll(building_t const &building, point &pos, point const &p_last, float radius,
		vect_room_object_t::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix) const;
	bool check_sphere_coll_walls_elevators_doors(building_t const &building, point &pos, point const &p_last, float radius,
		float wall_test_extra_z, bool check_open_doors, vector3d *cnorm) const;
	bool line_coll(building_t const &building, point const &p1, point const &p2, point &p_int) const;
	point find_closest_pt_on_obj_to_pos(building_t const &building, point const &pos, float pad_dist, bool no_ceil_floor) const;
	void update_dynamic_draw_data() {assert(room_geom); room_geom->update_dynamic_draw_data();}
	void get_avoid_cubes(vect_cube_t &avoid, float z1, float z2, float floor_thickness, bool same_as_player) const;
};

struct building_stats_t {
	unsigned nbuildings, nparts, ndetails, ntquads, ndoors, ninterior, nrooms, nceils, nfloors, nwalls, nrgeom, nobjs, nverts;
	building_stats_t() : nbuildings(0), nparts(0), ndetails(0), ntquads(0), ndoors(0), ninterior(0), nrooms(0), nceils(0), nfloors(0), nwalls(0), nrgeom(0), nobjs(0), nverts(0) {}
};

struct building_loc_t {
	int part_ix, room_ix, stairs_ix; // -1 is not contained; what about elevator_ix?
	unsigned floor_ix; // global for this building, rather than the current part/room
	building_loc_t() : part_ix(-1), room_ix(-1), stairs_ix(-1), floor_ix(0) {}
	bool operator==(building_loc_t const &loc) const {return (same_room_floor(loc) && part_ix == loc.part_ix && stairs_ix == loc.stairs_ix);}
	bool same_room_floor(building_loc_t const &loc) const {return (room_ix == loc.room_ix && floor_ix == loc.floor_ix);}
};

struct building_dest_t : public building_loc_t {
	int building_ix;
	point pos;
	building_dest_t() : building_ix(-1) {}
	building_dest_t(building_loc_t const &b, point const &pos_, int bix=-1) : building_loc_t(b), building_ix(bix), pos(pos_) {}
	bool is_valid() const {return (building_ix >= 0 && part_ix >= 0 && room_ix >= 0);}
};

enum {AI_STOP=0, AI_WAITING, AI_NEXT_PT, AI_BEGIN_PATH, AI_AT_DEST, AI_MOVING};
enum {GOAL_TYPE_NONE=0, GOAL_TYPE_ROOM, GOAL_TYPE_PLAYER, GOAL_TYPE_SOUND};

struct building_ai_state_t {
	bool is_first_path, on_new_path_seg;
	uint8_t goal_type;
	int cur_room, dest_room; // Note: -1 is unassigned
	vector<point> path; // stored backwards, next point on path is path.back()

	building_ai_state_t() : is_first_path(1), on_new_path_seg(0), goal_type(GOAL_TYPE_NONE), cur_room(-1), dest_room(-1) {}
	void next_path_pt(pedestrian_t &person, bool same_floor, bool starting_path);
};

struct colored_cube_t;
typedef vector<colored_cube_t> vect_colored_cube_t;
class cube_bvh_t;
class building_indir_light_mgr_t;

struct building_t : public building_geom_t {

	unsigned mat_ix;
	uint8_t hallway_dim, real_num_parts, roof_type; // main hallway dim: 0=x, 1=y, 2=none
	uint8_t street_dir; // encoded as 2*dim + dir + 1; 0 is unassigned
	int8_t open_door_ix, basement_part_ix, has_chimney; // has_chimney: 0=none, 1=interior, 2=exterior with fireplace
	bool is_house, has_garage, has_shed, has_int_garage, has_courtyard, has_complex_floorplan, has_helipad, has_ac, has_int_fplace;
	colorRGBA side_color, roof_color, detail_color, door_color, wall_color;
	cube_t bcube, pri_hall, driveway, porch, assigned_plot;
	vect_cube_t parts, fences;
	vect_roof_obj_t details; // cubes on the roof - antennas, AC units, etc.
	vector<tquad_with_ix_t> roof_tquads, doors;
	std::shared_ptr<building_interior_t> interior;
	vertex_range_t ext_side_qv_range;
	point tree_pos; // (0,0,0) is unplaced/no tree
	float ao_bcz2, ground_floor_z1;

	friend class building_indir_light_mgr_t;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_), hallway_dim(2), real_num_parts(0), roof_type(ROOF_TYPE_FLAT), street_dir(0), open_door_ix(-1),
		basement_part_ix(-1), has_chimney(0), is_house(0), has_garage(0), has_shed(0), has_int_garage(0), has_courtyard(0), has_complex_floorplan(0),
		has_helipad(0), has_ac(0), has_int_fplace(0), side_color(WHITE), roof_color(WHITE), detail_color(BLACK), door_color(WHITE), wall_color(WHITE),
		ao_bcz2(0.0), ground_floor_z1(0.0) {}
	building_t(building_geom_t const &bg) : building_geom_t(bg), mat_ix(0), hallway_dim(2), real_num_parts(0), roof_type(ROOF_TYPE_FLAT), street_dir(0), open_door_ix(-1),
		basement_part_ix(-1), has_chimney(0), is_house(0), has_garage(0), has_shed(0), has_int_garage(0), has_courtyard(0), has_complex_floorplan(0), has_helipad(0),
		has_ac(0), has_int_fplace(0), side_color(WHITE), roof_color(WHITE), detail_color(BLACK), door_color(WHITE), wall_color(WHITE), ao_bcz2(0.0), ground_floor_z1(0.0) {}
	static float get_scaled_player_radius();
	static float get_min_front_clearance() {return 2.05f*get_scaled_player_radius();} // slightly larger than the player diameter
	bool is_valid() const {return !bcube.is_all_zeros();}
	bool has_interior () const {return bool(interior);}
	bool has_room_geom() const {return (has_interior() && interior->room_geom);}
	bool has_sec_bldg () const {return (has_garage || has_shed);}
	bool has_pri_hall () const {return (hallway_dim <= 1);} // otherswise == 2
	bool has_basement () const {return (basement_part_ix >= 0);}
	bool has_driveway () const {return !driveway.is_all_zeros();}
	bool has_a_garage () const {return (has_garage || has_int_garage);} // external or internal
	bool enable_driveway_coll() const {return !is_rotated();} // no collision with rotated driveways/porches for now
	bool is_basement(vect_cube_t::const_iterator it) const {return (int(it - parts.begin()) == basement_part_ix);}
	bool is_pos_in_basement(point const &pos) const {return (has_basement() && parts[basement_part_ix].contains_pt(pos));};
	colorRGBA get_avg_side_color  () const {return side_color  .modulate_with(get_material().side_tex.get_avg_color());}
	colorRGBA get_avg_roof_color  () const {return roof_color  .modulate_with(get_material().roof_tex.get_avg_color());}
	colorRGBA get_avg_detail_color() const {return detail_color.modulate_with(get_material().roof_tex.get_avg_color());}
	building_mat_t const &get_material() const;
	bool has_windows() const {return get_material().add_windows;}
	float get_floor_thick_val() const {return (is_house ? FLOOR_THICK_VAL_HOUSE : (has_windows() ? FLOOR_THICK_VAL_OFFICE : FLOOR_THICK_VAL_WINDOWLESS));}
	float get_elevator_fc_thick_scale() const {return 1.005*0.5*get_floor_thick_val();}
	float get_window_vspace  () const {return get_material().get_floor_spacing();}
	float get_floor_thickness() const {return get_floor_thick_val()*get_window_vspace();}
	float get_fc_thickness   () const {return 0.5*get_floor_thickness();} // floor/ceiling thickness
	float get_wall_thickness () const {return WALL_THICK_VAL*get_window_vspace();}
	float get_trim_thickness () const {return 0.1*get_wall_thickness();}
	float get_trim_height    () const {return 0.04*get_window_vspace();}
	float get_door_height    () const {return 0.95f*(get_window_vspace() - get_floor_thickness());} // set height based on window spacing, 95% of ceiling height (may be too large)
	float get_doorway_width  () const;
	float get_ground_floor_z_thresh() const {return (ground_floor_z1 + 0.25f*get_window_vspace());} // for rats
	unsigned get_person_capacity_mult() const;
	void gen_rotation(rand_gen_t &rgen);
	void maybe_inv_rotate_point(point &p) const {if (is_rotated()) {do_xy_rotate_inv(bcube.get_cube_center(), p);}} // inverse rotate - negate the sine term
	void maybe_inv_rotate_pos_dir(point &pos, vector3d &dir) const;
	void set_z_range(float z1, float z2);
	bool check_part_contains_pt_xy(cube_t const &part, point const &pt, vector<point> &points) const;
	bool check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const;
	vect_cube_t::const_iterator get_real_parts_end() const {return (parts.begin() + real_num_parts);}
	vect_cube_t::const_iterator get_real_parts_end_inc_sec() const {return (get_real_parts_end() + has_sec_bldg());}
	cube_t const &get_sec_bldg () const {assert(has_sec_bldg()); assert(real_num_parts < parts.size()); return parts[real_num_parts];}
	cube_t const &get_chimney  () const {assert(has_chimney == 2 && parts.size() > 1); return parts.back();}
	cube_t const &get_fireplace() const {assert(has_chimney == 2 && parts.size() > 2); return parts[parts.size()-2];}
	void end_add_parts() {assert(parts.size() < 256); real_num_parts = uint8_t(parts.size());}
	cube_t get_coll_bcube() const;

	bool check_sphere_coll(point const &pos, float radius, bool xy_only, vector<point> &points, vector3d *cnorm=nullptr) const {
		point pos2(pos);
		return check_sphere_coll(pos2, pos, vect_cube_t(), zero_vector, radius, xy_only, points, cnorm);
	}
	bool check_sphere_coll(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, vector3d const &xlate, float radius, bool xy_only,
		vector<point> &points, vector3d *cnorm=nullptr, bool check_interior=0) const;
	bool check_sphere_coll_interior(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, float radius, bool xy_only, vector3d *cnorm=nullptr) const;
	unsigned check_line_coll(point const &p1, point const &p2, float &t, vector<point> &points, bool occlusion_only=0, bool ret_any_pt=0, bool no_coll_pt=0) const;
	bool check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const;
	bool check_point_xy_in_part(point const &pos) const;
	bool ray_cast_exterior_walls(point const &p1, point const &p2, vector3d &cnorm, float &t) const;
	bool ray_cast_interior(point const &pos, vector3d const &dir, cube_bvh_t const &bvh, point &cpos, vector3d &cnorm, colorRGBA &ccolor) const;
	void create_building_volume_light_texture(unsigned bix, point const &target, unsigned &tid) const;
	bool ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const;
	void calc_bcube_from_parts();
	void adjust_part_zvals_for_floor_spacing(cube_t &c) const;
	void gen_geometry(int rseed1, int rseed2);
	cube_t place_door(cube_t const &base, bool dim, bool dir, float door_height, float door_center, float door_pos,
		float door_center_shift, float width_scale, bool can_fail, bool opens_up, rand_gen_t &rgen) const;
	void gen_house(cube_t const &base, rand_gen_t &rgen);
	bool maybe_add_house_driveway(cube_t const &plot, cube_t &ret, unsigned building_ix) const;
	bool get_power_point(vector<point> &ppts) const;
	void add_solar_panels(rand_gen_t &rgen);
	bool add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_office_building, bool roof_access=0);
	float gen_peaked_roof(cube_t const &top_, float peak_height, bool dim, float extend_to, float max_dz, unsigned skip_side_tri);
	float gen_hipped_roof(cube_t const &top_, float peak_height, float extend_to);
	void place_roof_ac_units(unsigned num, float sz_scale, cube_t const &bounds, vect_cube_t const &avoid, bool avoid_center, rand_gen_t &rgen);
	void add_roof_walls(cube_t const &c, float wall_width, bool overlap_corners, cube_t out[4]);
	void gen_details(rand_gen_t &rgen, bool is_rectangle);
	cube_t get_helipad_bcube() const;
	int get_num_windows_on_side(float xy1, float xy2) const;
	float get_window_h_border() const;
	float get_window_v_border() const;
	float get_hspacing_for_part(cube_t const &part, bool dim) const;
	bool interior_enabled() const;
	void gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes);
	int maybe_assign_interior_garage(bool &gdim, bool &gdir);
	void add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned part_ix, unsigned num_floors,
		unsigned rooms_start, bool use_hallway, bool first_part_this_stack, float window_hspacing[2], float window_border);
	void connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part);
	void gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes, unsigned building_ix);
	void add_stairs_and_elevators(rand_gen_t &rgen);
	int get_ext_door_dir(cube_t const &door_bcube, bool dim) const;
	void add_sign_by_door(tquad_with_ix_t const &door, bool outside, std::string const &text, colorRGBA const &color, bool emissive);
	void add_doorbell(tquad_with_ix_t const &door);
	void add_exterior_door_signs(rand_gen_t &rgen);
	void gen_building_doors_if_needed(rand_gen_t &rgen);
	void maybe_add_special_roof(rand_gen_t &rgen);
	void gen_sloped_roof(rand_gen_t &rgen, cube_t const &top);
	void add_roof_to_bcube();
	void gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax);
	void get_all_drawn_verts(building_draw_t &bdraw, bool get_exterior, bool get_interior, bool get_int_ext_walls);
	void get_all_drawn_window_verts(building_draw_t &bdraw, bool lights_pass=0, float offset_scale=1.0, point const *const only_cont_pt_in=nullptr) const;
	void get_all_drawn_window_verts_as_quads(vect_vnctcc_t &verts) const;
	bool get_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, point const &pos, float dist);
	void get_all_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, vector<point> const &pts, float dist);
	void player_not_near_building() {register_open_ext_door_state(-1);}
	int find_ext_door_close_to_point(tquad_with_ix_t &door, point const &pos, float dist) const;
	bool get_building_door_pos_closest_to(point const &target_pos, point &door_pos) const;
	void get_split_int_window_wall_verts(building_draw_t &bdraw_front, building_draw_t &bdraw_back, point const &only_cont_pt_in, bool make_all_front=0) const;
	void get_ext_wall_verts_no_sec(building_draw_t &bdraw) const;
	void write_basement_entrance_depth_pass(shader_t &s) const;
	void add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, int ped_ix, occlusion_checker_noncity_t &oc, vect_cube_t &ped_bcubes, cube_t &lights_bcube);
	bool toggle_room_light(point const &closest_to, bool sound_from_closest_to=0, int room_id=-1, bool inc_lamps=1, bool closet_light=0);
	void toggle_light_object(room_object_t const &light, point const &sound_pos);
	bool apply_player_action_key(point const &closest_to_in, vector3d const &in_dir_in, int mode, bool check_only=0);
	bool move_nearest_object(point const &at_pos, vector3d const &in_dir, float range, int mode);
	bool interact_with_object(unsigned obj_ix, point const &int_pos, point const &query_ray_end, vector3d const &int_dir);
	bool adjust_blinds_state(unsigned obj_ix);
	void add_box_contents(room_object_t const &box);
	void toggle_door_state(unsigned door_ix, bool player_in_this_building, bool by_player, float zval);
	bool set_room_light_state_to(room_t const &room, float zval, bool make_on);
	void set_obj_lit_state_to(unsigned room_id, float light_z2, bool lit_state);
	bool player_pickup_object(point const &at_pos, vector3d const &in_dir);
	void register_player_enter_building() const;
	void register_player_exit_building () const;
	bool check_for_wall_ceil_floor_int(point const &p1, point const &p2) const;
	bool maybe_use_last_pickup_room_object(point const &player_pos);
	bool maybe_update_tape(point const &player_pos, bool end_of_tape);
	void handle_vert_cylin_tape_collision(point &cur_pos, point const &prev_pos, float z1, float z2, float radius, bool is_player) const;
	void draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, occlusion_checker_noncity_t &oc, vector3d const &xlate,
		unsigned building_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building);
	void gen_and_draw_room_geom(brg_batch_draw_t *bbd, shader_t &s, occlusion_checker_noncity_t &oc, vector3d const &xlate, vect_cube_t &ped_bcubes,
		unsigned building_ix, int ped_ix, bool shadow_only, bool reflection_pass, unsigned inc_small, bool player_in_building);
	void add_split_roof_shadow_quads(building_draw_t &bdraw) const;
	void clear_room_geom(bool force);
	void update_grass_exclude_at_pos(point const &pos, vector3d const &xlate, bool camera_in_building) const;
	void update_stats(building_stats_t &s) const;
	bool are_rooms_connected_without_using_room(unsigned room1, unsigned room2, unsigned room_exclude) const;
	bool is_room_adjacent_to_ext_door(cube_t const &room, bool front_door_only=0) const;
	room_t const &get_room(unsigned room_ix) const {assert(interior && room_ix < interior->rooms.size()); return interior->rooms[room_ix];}
	point get_center_of_room(unsigned room_ix) const {return get_room(room_ix).get_cube_center();}

	// building AI people
	unsigned count_connected_room_components();
	bool place_person(point &ppos, float radius, rand_gen_t &rgen) const;
	int ai_room_update(building_ai_state_t &state, rand_gen_t &rgen, vector<pedestrian_t> &people, float delta_dir, unsigned person_ix, bool stay_on_one_floor);
private:
	void build_nav_graph() const;
	bool is_valid_ai_placement(point const &pos, float radius) const;
	bool choose_dest_goal(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const;
	int  choose_dest_room(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const;
	void get_avoid_cubes(float zval, float height, float radius, vect_cube_t &avoid, bool following_player) const;
	bool find_route_to_point(pedestrian_t const &person, float radius, bool is_first_path, bool is_moving_target, bool following_player, vector<point> &path) const;
	void find_nearest_stairs(point const &p1, point const &p2, vector<unsigned> &nearest_stairs, bool straight_only, int part_ix=-1) const;
	void ai_room_lights_update(building_ai_state_t const &state, pedestrian_t const &person, vector<pedestrian_t> const &people, unsigned person_ix);
	void move_person_to_not_collide(pedestrian_t &person, pedestrian_t const &other, point const &new_pos, float rsum, float coll_dist) const;

	// animals
public:
	void update_animals(point const &camera_bs, unsigned building_ix, int ped_ix);
	void get_objs_at_or_below_ground_floor(vect_room_object_t &ret) const;
private:
	point gen_rat_pos(float radius, rand_gen_t &rgen) const;
	bool add_rat(point const &pos, float hlength, vector3d const &dir, point const &placed_from);
	bool is_rat_inside_building(point const &pos, float xy_pad, float hheight) const;
	void update_rat(rat_t &rat, point const &camera_bs, int ped_ix, float timestep, float &max_xmove, bool can_attack_player, rand_gen_t &rgen) const;
	void scare_rat(rat_t &rat, point const &camera_bs, int ped_ix) const;
	void scare_rat_at_pos(rat_t &rat, point const &scare_pos, float amount, bool by_sight) const;
	bool check_line_coll_expand(point const &p1, point const &p2, float radius, float height) const;
	bool check_line_of_sight_large_objs(point const &p1, point const &p2) const;
	bool check_and_handle_dynamic_obj_coll(point &pos, float radius, float height, point const &camera_bs) const;
	bool get_begin_end_room_objs_on_ground_floor(float zval, vect_room_object_t::const_iterator &b, vect_room_object_t::const_iterator &e) const;

public:
	int get_room_containing_pt(point const &pt) const;
	bool room_containing_pt_has_stairs(point const &pt) const;
	bool maybe_teleport_to_screenshot() const;
	bool place_obj_along_wall(room_object type, room_t const &room, float height, vector3d const &sz_scale, rand_gen_t &rgen,
		float zval, unsigned room_id, float tot_light_amt, cube_t const &place_area, unsigned objs_start, float front_clearance=0.0,
		unsigned pref_orient=4, bool pref_centered=0, colorRGBA const &color=WHITE, bool not_at_window=0, room_obj_shape shape=SHAPE_CUBE);
	bool place_model_along_wall(unsigned model_id, room_object type, room_t const &room, float height, rand_gen_t &rgen,
		float zval, unsigned room_id, float tot_light_amt, cube_t const &place_area, unsigned objs_start, float front_clearance=0.0,
		unsigned pref_orient=4, bool pref_centered=0, colorRGBA const &color=WHITE, bool not_at_window=0);
	int check_valid_picture_placement(room_t const &room, cube_t const &c, float width, float zval, bool dim, bool dir, unsigned objs_start) const;
	void update_player_interact_objects(point const &player_pos, int first_ped_ix);
	bool line_intersect_walls(point const &p1, point const &p2) const;
	bool is_obj_pos_valid(room_object_t const &obj, bool keep_in_room, bool allow_block_door, bool check_stairs) const;
	bool is_rot_cube_visible(cube_t const &c, vector3d const &xlate) const;
	bool is_cube_face_visible_from_pt(cube_t const &c, point const &p, unsigned dim, bool dir) const;
	bool check_obj_occluded(cube_t const &c, point const &viewer, occlusion_checker_noncity_t &oc, bool reflection_pass, bool c_is_building_part=0) const;
	bool is_entire_building_occluded(point const &viewer, occlusion_checker_noncity_t &oc) const;
	template<typename T> void add_door_verts(cube_t const &D, T &drawer, uint8_t door_type,
		bool dim, bool dir, bool opened, bool opens_out, bool exterior, bool on_stairs=0, bool hinge_side=0) const;
	tquad_with_ix_t set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj,
		bool exterior, bool opened, bool opens_out, bool opens_up, bool swap_sides) const;
	tquad_with_ix_t set_interior_door_from_cube(door_t const &door) const;
	void invalidate_nav_graph();
	point local_to_camera_space(point const &pos) const;
	void play_door_open_close_sound(point const &pos, bool open, float gain=1.0, float pitch=1.0) const;
	void maybe_gen_chimney_smoke() const;
	cube_t get_part_containing_pt(point const &pt) const;
	void print_building_manifest() const;
private:
	void finish_gen_geometry(rand_gen_t &rgen, bool has_overlapping_cubes);
	bool add_outdoor_ac_unit(rand_gen_t &rgen);
	bool add_chimney(cube_t const &part, bool dim, bool dir, float chimney_dz, int garage_room, rand_gen_t &rgen);
	void gen_interior_int(rand_gen_t &rgen, bool has_overlapping_cubes);
	void maybe_add_basement(rand_gen_t rgen);
	void clip_cube_to_parts(cube_t &c, vect_cube_t &cubes) const;
	bool move_sphere_to_valid_part(point &pos, point const &p_last, float radius) const;
	cube_t get_walkable_room_bounds(room_t const &room) const;
	bool is_cube_contained_in_parts(cube_t const &c) const;
	void expand_ground_floor_cube(cube_t &cube, cube_t const &skip=cube_t()) const;
	void get_exclude_cube(point const &pos, cube_t const &skip, cube_t &exclude, bool camera_in_building) const;
	void move_door_to_other_side_of_wall(tquad_with_ix_t &door, float dist_mult, bool invert_normal) const;
	void clip_door_to_interior(tquad_with_ix_t &door, bool clip_to_floor) const;
	void cut_holes_for_ext_doors(building_draw_t &bdraw, point const &contain_pt, unsigned draw_parts_mask) const;
	bool is_valid_door_pos(cube_t const &door, float door_width, bool dim) const;
	bool is_cube_close_to_exterior_doorway(cube_t const &c, float dmin=0.0, bool inc_open=0) const;
	bool is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin=0.0, bool inc_open=0, bool check_open_dir=0) const;
	bool is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, float room_pad=0.0) const;
	bool check_cube_intersect_walls(cube_t const &c) const;
	bool check_cube_contained_in_part(cube_t const &c) const;
	bool is_valid_stairs_elevator_placement(cube_t const &c, float pad, bool check_walls=1) const;
	bool clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const;
	unsigned add_room(cube_t const &room, unsigned part_id, unsigned num_lights, bool is_hallway, bool is_office, bool is_sec_bldg=0);
	void add_or_extend_elevator(elevator_t const &elevator, bool add);
	void remove_intersecting_roof_cubes(cube_t const &c);
	bool overlaps_other_room_obj(cube_t const &c, unsigned objs_start=0, bool check_all=0) const;
	bool overlaps_any_placed_obj(cube_t const &c) const;
	int classify_room_wall(room_t const &room, float zval, bool dim, bool dir, bool ret_sep_if_part_int_part_ext) const;
	unsigned count_ext_walls_for_room(room_t const &room, float zval) const;
	bool room_has_stairs_or_elevator(room_t const &room, float zval, unsigned floor) const;
	bool is_room_office_bathroom(room_t const &room, float zval, unsigned floor) const;
	int gather_room_placement_blockers(cube_t const &room, unsigned objs_start, vect_cube_t &blockers, bool inc_open_doors=1, bool ignore_chairs=0) const;
	vect_door_stack_t &get_doorways_for_room(room_t const &room, float zval) const;
	bool add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos,
		colorRGBA const &chair_color, bool dim, bool dir, float tot_light_amt, bool office_chair_model);
	unsigned add_table_and_chairs(rand_gen_t rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
		point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt);
	void shorten_chairs_in_region(cube_t const &region, unsigned objs_start);
	void add_trashcan_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool check_last_obj);
	bool add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	bool add_desk_to_room    (rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, colorRGBA const &chair_color,
		float zval, unsigned room_id, unsigned floor, float tot_light_amt, bool is_basement);
	bool create_office_cubicles(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	bool check_valid_closet_placement(cube_t const &c, room_t const &room, unsigned objs_start, unsigned bed_ix, float min_bed_space=0.0) const;
	bool add_bedroom_objs    (rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, unsigned floor,
		float tot_light_amt, unsigned objs_start, bool room_is_lit, bool is_basement, light_ix_assign_t &light_ix_assign);
	bool add_bed_to_room     (rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, float tot_light_amt, unsigned floor);
	bool maybe_add_fireplace_to_room(room_t const &room, vect_cube_t &blockers, float zval, unsigned room_id, float tot_light_amt);
	float add_flooring       (room_t const &room, float &zval, unsigned room_id, float tot_light_amt);
	bool add_bathroom_objs   (rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt,
		unsigned objs_start, unsigned floor, bool is_basement, unsigned &added_bathroom_objs_mask);
	bool add_tp_roll         (cube_t const &room, unsigned room_id, float tot_light_amt, bool dim, bool dir, float length, float zval, float wall_pos, bool check_valid_pos=0);
	bool divide_bathroom_into_stalls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned floor);
	bool add_kitchen_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool allow_adj_ext_door);
	bool add_livingroom_objs (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	void add_diningroom_objs (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_library_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	bool add_storage_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	bool add_basement_utility_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool add_laundry_objs    (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned &added_bathroom_objs_mask);
	void add_pri_hall_objs   (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt);
	void add_parking_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start);
	void place_book_on_obj   (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, bool use_dim_dir);
	bool place_bottle_on_obj (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid);
	bool place_plant_on_obj  (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid);
	bool place_laptop_on_obj (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid, bool use_dim_dir);
	bool place_plate_on_obj  (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid);
	bool place_cup_on_obj    (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid);
	bool add_rug_to_room     (rand_gen_t rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start);
	bool hang_pictures_in_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement);
	void add_plants_to_room  (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num);
	void add_boxes_to_room   (rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned max_num);
	void add_light_switches_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor);
	void add_outlets_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor, bool is_basement);
	bool place_eating_items_on_table(rand_gen_t &rgen, unsigned table_obj_id);
	void place_objects_onto_surfaces(rand_gen_t rgen, room_t const &room, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned floor, bool is_basement);
	bool can_be_bedroom_or_bathroom(room_t const &room, unsigned floor, bool skip_conn_check=0) const;
	bool can_be_bathroom(room_t const &room) const;
	bool find_mirror_in_room(unsigned room_id, vector3d const &xlate, bool check_visibility) const;
	bool find_mirror_needing_reflection(vector3d const &xlate) const;
	tquad_with_ix_t const &find_main_roof_tquad(rand_gen_t &rgen, bool skip_if_has_other_obj) const;
	void maybe_add_fire_escape(rand_gen_t &rgen);
	void add_extra_obj_slots();
	void add_wall_and_door_trim_if_needed();
	void add_wall_and_door_trim();
	unsigned count_num_int_doors(room_t const &room) const;
	bool check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const;
	void split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen);
	bool test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const;
	void gather_interior_cubes(vect_colored_cube_t &cc) const;
	void order_lights_by_priority(point const &target, vector<unsigned> &light_ids) const;
	bool is_light_occluded(point const &lpos, point const &camera_bs) const;
	void clip_ray_to_walls(point const &p1, point &p2) const;
	void refine_light_bcube(point const &lpos, float light_radius, cube_t const &room, cube_t &light_bcube) const;
	cube_t get_rotated_bcube(cube_t const &c) const;
	cube_t get_part_for_room(room_t const &room) const {assert(room.part_id < parts.size()); return parts[room.part_id];}
	bool are_parts_stacked(cube_t const &p1, cube_t const &p2) const;
	void add_window_coverings(cube_t const &window, bool dim, bool dir);
	void add_window_blinds(cube_t const &window, bool dim, bool dir, unsigned room_ix, unsigned floor);
	void add_bathroom_window(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor);
	int get_room_id_for_window(cube_t const &window, bool dim, bool dir, bool &is_split) const;
	void register_open_ext_door_state(int door_ix);
	void add_interior_door(door_t &door, bool is_bathroom=0);
	void add_interior_door_for_floor(door_t &door, bool is_bathroom);
	void remove_section_from_cube_and_add_door(cube_t &c, cube_t &c2, float v1, float v2, bool xy, bool open_dir, bool is_bathroom=0);
	void insert_door_in_wall_and_add_seg(cube_t &wall, float v1, float v2, bool dim, bool open_dir, bool keep_high_side, bool is_bathroom=0);
	unsigned get_floor_for_zval(float zval) const {return unsigned((zval - bcube.z1())/get_window_vspace());}
	building_loc_t get_building_loc_for_pt(point const &pt) const;
	void register_player_in_building(point const &camera_bs, unsigned building_id) const;
	bool same_room_and_floor_as_player(building_ai_state_t const &state, pedestrian_t const &person) const;
	bool is_player_visible(building_ai_state_t const &state, pedestrian_t const &person, unsigned vis_test) const;
	bool can_target_player(building_ai_state_t const &state, pedestrian_t const &person) const;
	bool need_to_update_ai_path(building_ai_state_t const &state, pedestrian_t const &person) const;
	void set_bcube_from_rotated_cube(cube_t const &bc);
	bool apply_paint(point const &pos, vector3d const &dir, colorRGBA const &color, room_object const obj_type) const;
	bool apply_toilet_paper(point const &pos, vector3d const &dir, float half_width);
	void register_button_event(room_object_t const &button);
	bool get_zval_of_floor(point const &pos, float radius, float &zval) const;
	bool get_zval_for_obj_placement(point const &pos, float radius, float &zval, bool add_z_bias) const;
	void add_blood_decal(point const &pos);
	void play_tape_sound(point const &sound_pos, float sound_gain) const;
public:
	// ray queries
	bool check_line_intersect_doors(point const &p1, point const &p2, bool inc_open=0) const;
	bool is_pt_visible(point const &p1, point const &p2) const;
	bool is_sphere_visible(point const &center, float radius, point const &pt) const;
	bool is_pt_lit(point const &pt) const;
	bool is_sphere_lit(point const &center, float radius) const;
};

struct vect_building_t : public vector<building_t> {
	vect_cube_with_ix_t building_to_person;
	void ai_room_update(vector<building_ai_state_t> &ai_state, vector<pedestrian_t> &people, unsigned p_start, unsigned p_end, float delta_dir, rand_gen_t &rgen);
};

struct building_draw_utils {
	static void calc_normals(building_geom_t const &bg, vector<vector3d> &nv, unsigned ndiv);
	static void calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, vector<point> &pts, float expand=0.0);
};

class city_lights_manager_t {
protected:
	cube_t lights_bcube;
	float light_radius_scale, dlight_add_thresh;
	bool prev_had_lights;
public:
	city_lights_manager_t() : lights_bcube(all_zeros), light_radius_scale(1.0), dlight_add_thresh(0.0), prev_had_lights(0) {}
	virtual ~city_lights_manager_t() {}
	cube_t get_lights_bcube() const {return lights_bcube;}
	void add_player_flashlight(float radius_scale);
	void tighten_light_bcube_bounds(vector<light_source> const &lights);
	void clamp_to_max_lights(vector3d const &xlate, vector<light_source> &lights);
	bool begin_lights_setup(vector3d const &xlate, float light_radius, vector<light_source> &lights);
	void finalize_lights(vector<light_source> &lights);
	void setup_shadow_maps(vector<light_source> &light_sources, point const &cpos, unsigned max_smaps);
	virtual bool enable_lights() const = 0;
};

struct building_place_t {
	point p;
	unsigned bix;
	building_place_t() : bix(0) {}
	building_place_t(point const &p_, unsigned bix_) : p(p_), bix(bix_) {}
	bool operator<(building_place_t const &p) const {return (bix < p.bix);} // only compare building index
};
typedef vector<building_place_t> vect_building_place_t;

struct ped_draw_vars_t {
	building_t const &building;
	occlusion_checker_noncity_t &oc;
	shader_t &s;
	vector3d const &xlate;
	unsigned bix;
	bool shadow_only, reflection_pass;

	ped_draw_vars_t(building_t const &b, occlusion_checker_noncity_t &oc_, shader_t &s_, vector3d const &x, unsigned bix_, bool so, bool rp)
		: building(b), oc(oc_), s(s_), xlate(x), bix(bix_), shadow_only(so), reflection_pass(rp) {}
};

class water_sound_manager_t {
	point const camera_bs;
	point closest; // in camera space
	float dmin_sq;
public:
	water_sound_manager_t(point const &camera_bs_) : camera_bs(camera_bs_), dmin_sq(0.0) {}
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

template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &bcubes) { // T must derive from cube_t
	for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {if (c->intersects(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_no_adj(cube_t const &bcube, vector<T> const &bcubes) { // T must derive from cube_t
	for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {if (c->intersects_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> bool has_bcube_int_xy_no_adj(cube_t const &bcube, vector<T> const &bcubes) { // T must derive from cube_t
	for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {if (c->intersects_xy_no_adj(bcube)) return 1;}
	return 0;
}
template<typename T> cube_t get_cube_height_radius(point const &center, T radius, float height) { // T can be float or vector3d
	cube_t c(center);
	c.expand_by_xy(radius);
	c.z2() += height;
	return c;
}

inline point get_camera_building_space() {return (get_camera_pos() - get_tiled_terrain_model_xlate());}
inline void set_cube_zvals(cube_t &c, float z1, float z2) {c.z1() = z1; c.z2() = z2;}
inline float get_tc_leg_width(cube_t const &c, float width) {return 0.5f*width*(c.dx() + c.dy());} // make legs square

void get_city_plot_zones(vect_city_zone_t &zones);
void get_city_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_city_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state);
bool city_single_cube_visible_check(point const &pos, cube_t const &c);
cube_t get_building_lights_bcube();
unsigned get_street_dir(cube_t const &inner, cube_t const &outer);
cube_t get_open_closet_door(room_object_t const &c, cube_t const &closed_door);
float get_closet_wall_thickness(room_object_t const &c);
void get_closet_cubes(room_object_t const &c, cube_t cubes[5], bool for_collision=0);
void get_bed_cubes   (room_object_t const &c, cube_t cubes[6]);
void get_table_cubes (room_object_t const &c, cube_t cubes[5]);
void get_chair_cubes (room_object_t const &c, cube_t cubes[3]);
void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]);
void get_bookcase_cubes(room_object_t const &c, cube_t &top, cube_t &middle, cube_t &back, cube_t lr[2], bool no_shelves=0, float sides_scale=1.0);
float get_drawer_cubes(room_object_t const &c, vect_cube_t &drawers, bool front_only, bool inside_only);
float get_cabinet_doors(room_object_t const &c, vect_cube_t &doors);
void get_cabinet_or_counter_doors(room_object_t const &c, vect_cube_t &doors);
bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher);
room_object_t split_cabinet_at_dishwasher(room_object_t &cabinet, cube_t const &dishwasher);
room_object_t get_dresser_middle(room_object_t const &c);
room_object_t get_desk_drawers_part(room_object_t const &c);
cube_t get_elevator_car_panel(room_object_t const &c, float fc_thick_scale);
void set_rand_pos_for_sz(cube_t &c, bool dim, float length, float width, rand_gen_t &rgen);
template<typename T> bool has_bcube_int_xy(cube_t const &bcube, vector<T> const &bcubes, float pad_dist=0.0);
bool door_opens_inward(door_base_t const &door, cube_t const &room);
bool is_cube_close_to_door(cube_t const &c, float dmin, bool inc_open, cube_t const &door, unsigned check_dirs=2, unsigned open_dirs=2, bool allow_block_door=0);
void add_building_interior_lights(point const &xlate, cube_t &lights_bcube);
unsigned calc_num_floors(cube_t const &c, float window_vspacing, float floor_thickness);
void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim);
bool is_val_inside_window(cube_t const &c, bool dim, float val, float window_spacing, float window_border);
void subtract_cube_from_cube(cube_t const &c, cube_t const &s, vect_cube_t &out);
void subtract_cube_from_cube_inplace(cube_t const &s, vect_cube_t &cubes, unsigned &ix, unsigned &iter_end);
template<typename T> void subtract_cubes_from_cube(cube_t const &c, T const &sub, vect_cube_t &out, vect_cube_t &out2);
bool subtract_cube_from_cubes(cube_t const &s, vect_cube_t &cubes, vect_cube_t *holes=nullptr, bool clip_in_z=0, bool include_adj=0);
int get_rect_panel_tid();
int get_bath_wind_tid ();
int get_int_door_tid  ();
int get_normal_map_for_bldg_tid(int tid);
unsigned register_sign_text(std::string const &text);
void setup_building_draw_shader(shader_t &s, float min_alpha, bool enable_indir, bool force_tsl, bool use_texgen);
void rotate_verts(vector<rgeom_mat_t::vertex_t> &verts, building_t const &building);
void add_tquad_to_verts(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex,
	colorRGBA const &color, vect_vnctcc_t &verts, bool invert_tc_x=0, bool exclude_frame=0, bool no_tc=0, bool no_rotate=0);
void get_driveway_sphere_coll_cubes(point const &pos, float radius, bool xy_only, vect_cube_t &out);
bool have_buildings_ext_paint();
void draw_buildings_ext_paint();
void subtract_cube_xy(cube_t const &c, cube_t const &r, cube_t *out);
bool have_secondary_buildings();
bool get_building_door_pos_closest_to(unsigned building_id, point const &target_pos, point &door_pos);
void register_achievement(std::string const &str);
// functions in building_room_obj_expand.cc
point gen_xy_pos_in_area(cube_t const &S, vector3d const &sz, rand_gen_t &rgen);
point gen_xy_pos_in_area(cube_t const &S, float radius, rand_gen_t &rgen);
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen, bool place_at_z1=0);
// functions in building_interact.cc and building_gameplay.cc
void gen_sound_thread_safe(unsigned id, point const &pos, float gain=1.0, float pitch=1.0, float gain_scale=1.0, bool skip_if_already_playing=0);
inline void gen_sound_thread_safe_at_player(unsigned id, float gain=1.0, float pitch=1.0) {gen_sound_thread_safe(id, get_camera_pos(), gain, pitch);}
void register_building_sound(point const &pos, float volume);
void register_building_sound_at_player(float volume);
// functions in city_gen.cc
void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, int use_smap, int use_bmap,
	float min_alpha=0.0, bool force_tsl=0, float pcf_scale=1.0, bool use_texgen=0, bool indir_lighting=0, bool is_outside=1);
void enable_animations_for_shader(shader_t &s);
void setup_city_lights(vector3d const &xlate);
void draw_peds_in_building(int first_ped_ix, ped_draw_vars_t const &pdv); // from city_gen.cpp
void get_locations_of_peds_in_building(int first_ped_ix, vector<point> &locs); // from city_gen.cpp
void get_ped_bcubes_for_building(int first_ped_ix, vect_cube_t &bcubes, bool moving_only=0); // from city_gen.cpp
void register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity);
void draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only);
vector3d get_nom_car_size();
bool car_can_fit(cube_t const &c);
void draw_cars_in_garages(vector3d const &xlate, bool shadow_only);
void create_mirror_reflection_if_needed();
void draw_city_roads(int trans_op_mask, vector3d const &xlate);
void get_closest_dim_dir_xy(cube_t const &inner, cube_t const &outer, bool &dim, bool &dir);
bool check_city_tline_cube_intersect_xy(cube_t const &c);

