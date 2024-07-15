// 3D World - City Object Placement Header
// by Frank Gennari
// 08/14/21

#pragma once

#include "city.h"

void disable_hemi_lighting_pre_post(draw_state_t &dstate, bool shadow_only, bool is_post);

struct textured_mat_t {
	bool has_alpha_mask=0;
	int tid=-1, nm_tid=-1;
	colorRGBA color, map_color;
	string tex_name, nm_tex_name;

	textured_mat_t(string const &tn, string const &nm_tn, bool ham, colorRGBA const &c, colorRGBA const &mc) :
		has_alpha_mask(ham), color(c), map_color(mc), tex_name(tn), nm_tex_name(nm_tn) {}
	colorRGBA get_avg_color() const {return ((tid >= 0) ? texture_color(tid) : map_color).modulate_with(color);} // for overhead map
	void pre_draw (bool shadow_only);
	void post_draw(bool shadow_only);
};
struct plot_divider_type_t : public textured_mat_t {
	bool is_occluder;
	float wscale, hscale, tscale; // width, height, and texture scales
	plot_divider_type_t(string const &tn, string const &nm_tn, float ws, float hs, float ts, bool ic, bool ham, colorRGBA const &c, colorRGBA const &mc) :
		textured_mat_t(tn, nm_tn, ham, c, mc), is_occluder(ic), wscale(ws), hscale(hs), tscale(ts) {}
};
enum {DIV_WALL=0, DIV_FENCE, DIV_HEDGE, DIV_CHAINLINK, DIV_NUM_TYPES}; // types of plot dividers, with end terminator
enum {POOL_DECK_WOOD=0, POOL_DECK_CONCRETE, NUM_POOL_DECK_TYPES};
unsigned const NUM_POOL_DECK_PASSES(NUM_POOL_DECK_TYPES + 2); // {deck types}, roof, pillars

struct city_draw_qbds_t {
	quad_batch_draw qbd, untex_qbd, untex_spec_qbd, emissive_qbd;
	bool empty() const {return (qbd.empty() && untex_qbd.empty() && untex_spec_qbd.empty() && emissive_qbd.empty());}
	bool has_untex_verts() const {return (!untex_qbd.empty() || !untex_spec_qbd.empty());}
};

struct city_obj_t : public sphere_t {
	cube_t bcube;
	city_obj_t() {}
	city_obj_t(cube_t const &bcube_) : bcube(bcube_) {set_bsphere_from_bcube();}
	city_obj_t(point const &pos_, float radius_) : sphere_t(pos_, radius_) {}
	bool operator<(city_obj_t const &b) const {return (bcube.x1() < b.bcube.x1());} // sort by bcube x1
	cube_t const &get_outer_bcube() const {return bcube;}
	float get_bsphere_radius(bool shadow_only) const {return radius;}
	float get_overlay_radius() const {return radius;} // for overhead map mode
	void set_bsphere_from_bcube() {*((sphere_t *)this) = bcube.get_bsphere();}
	void set_bcube_from_vcylin(point const &base, float height, float xy_radius);
	static void pre_draw (draw_state_t &dstate, bool shadow_only) {} // nothing to do
	static void post_draw(draw_state_t &dstate, bool shadow_only) {} // nothing to do
	cube_t get_bird_bcube() const {return bcube;}
	bool check_point_contains_xy(point const &p) const {return bcube.contains_pt_xy(p);}
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const; // default bcube coll, can override in derived class
	void translate_dim(unsigned dim, float v) {pos[dim] += v; bcube.translate_dim(dim, v);}
};

struct oriented_city_obj_t : public city_obj_t {
	bool dim, dir;
	oriented_city_obj_t(bool dim_=0, bool dir_=0) : dim(dim_), dir(dir_) {}
	oriented_city_obj_t(cube_t const &bcube_, bool dim_=0, bool dir_=0) : city_obj_t(bcube_), dim(dim_), dir(dir_) {}
	oriented_city_obj_t(point const &pos_, float radius_, bool dim_=0, bool dir_=0) : city_obj_t(pos_, radius_), dim(dim_), dir(dir_) {}
	float get_length() const {return bcube.get_sz_dim( dim);}
	float get_depth () const {return bcube.get_sz_dim( dim);}
	float get_width () const {return bcube.get_sz_dim(!dim);}
	vector3d get_orient_dir() const;
};

struct model_city_obj_t : public oriented_city_obj_t {
	colorRGBA color=WHITE;
	float min_alpha=0.0;
	bool is_cylinder=0;

	model_city_obj_t(cube_t const &bcube_, bool dim_, bool dir_) : oriented_city_obj_t(bcube_, dim_, dir_) {}
	model_city_obj_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_id, bool is_cylinder_=0);
	virtual ~model_city_obj_t() {}
	virtual unsigned get_model_id() const = 0;
	virtual float get_xy_coll_radius() const {assert(is_cylinder); return 0.25*(bcube.dx() + bcube.dy());} // assume square-ish
	float get_overlay_radius() const {return model_city_obj_t::get_xy_coll_radius();} // for overhead map mode
	static void pre_draw (draw_state_t &dstate, bool shadow_only) {disable_hemi_lighting_pre_post(dstate, shadow_only, 0);}
	static void post_draw(draw_state_t &dstate, bool shadow_only) {disable_hemi_lighting_pre_post(dstate, shadow_only, 1);}
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only, animation_state_t *anim_state=nullptr) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct multi_model_city_obj_t : public model_city_obj_t {
	unsigned full_model_id=0;
	multi_model_city_obj_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_id, unsigned model_select, bool is_cylinder_=0);
	virtual unsigned get_model_id() const {return full_model_id;}
};

struct bench_t : public oriented_city_obj_t {
	bench_t(point const &pos_, float radius_, bool dim_, bool dir_);
	cube_t get_bird_bcube() const;
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct tree_planter_t : public city_obj_t {
	tree_planter_t(point const &pos_, float radius_, float height);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct trashcan_t : public city_obj_t {
	bool is_cylin;
	trashcan_t(point const &pos_, float radius_, float height, bool is_cylin_);
	float get_cylin_radius() const {assert(is_cylin); return 0.5*bcube.dx();}
	//cube_t get_bird_bcube() const {} // custom, not single cube; handled during object placement
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct fire_hydrant_t : public city_obj_t {
	float cylin_radius;
	vector3d orient;
	fire_hydrant_t(point const &pos_, float radius_, float height, vector3d const &orient_);
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct substation_t : public model_city_obj_t {
	substation_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id()) {}
	virtual unsigned get_model_id() const {return OBJ_MODEL_SUBSTATION;}
};

struct fountain_t : public multi_model_city_obj_t {
	fountain_t(point const &pos_, float height, unsigned model_select) :
		multi_model_city_obj_t(pos_, height, 0, 0, OBJ_MODEL_FOUNTAIN, model_select, 1) {} // dim=0, dir=0, is_cylinder=1
};

struct divider_t : public oriented_city_obj_t {
	unsigned type;
	uint8_t skip_dims;
	divider_t() : type(0), skip_dims(0) {}
	divider_t(cube_t const &c, unsigned type_, bool dim_, bool dir_, unsigned sd=0) :
		oriented_city_obj_t(c.get_cube_center(), c.get_bsphere_radius(), dim_, dir_), type(type_), skip_dims(sd) {bcube = c;}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct swimming_pool_t : public oriented_city_obj_t { // Note: dim and dir are used for the ladder and face toward the house
	colorRGBA color, wcolor; // wall color and water color
	bool above_ground, sloped;

	swimming_pool_t(cube_t const &c, colorRGBA const &color_, colorRGBA const &wcolor_, bool above_ground_, bool sloped_, bool dim_, bool dir_) :
		oriented_city_obj_t(c.get_cube_center(), c.get_bsphere_radius(), dim_, dir_), color(color_), wcolor(wcolor_), above_ground(above_ground_), sloped(sloped_) {bcube = c;}
	float get_radius() const {assert(above_ground); return 0.25f*(bcube.dx() + bcube.dy());}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct pool_ladder_t : public model_city_obj_t { // for in-ground pools
	pool_ladder_t(cube_t const &bcube_, bool dim_, bool dir_) : model_city_obj_t(bcube_, dim_, dir_) {}
	virtual unsigned get_model_id() const {return OBJ_MODEL_POOL_LAD;}
};

struct pool_deck_t : public oriented_city_obj_t {
	unsigned mat_id;
	cube_t base, roof;
	vect_cube_t pillars;

	pool_deck_t(cube_t const &base_, cube_t const &roof_, unsigned mat_id_, bool dim_, bool dir_);
	void calc_pillars(cube_t const &ladder);
	bool has_roof() const {return !roof.is_all_zeros();}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

class power_pole_t : public city_obj_t {
	struct wire_t {
		point pts[2], pole_base;
		wire_t(point const &p1, point const &p2) : pole_base(p1) {pts[0] = p1; pts[1] = p2;}
	};
	bool at_grid_edge, at_line_end[2], residential;
	uint8_t dims; // bit mask for direction the wires run
	float pole_radius, bsphere_radius, wires_offset, pole_spacing[2];
	point base, center; // base of the pole and center of wires/bcube
	cube_t bcube_with_wires;
	vector<wire_t> wires; // in addition to the normal connector wires

	float get_wire_radius    () const {return 0.08*pole_radius;}
	float get_bar_extend     () const {return 8.00*pole_radius;} // distance from the center that the wooden bar holding the wires extends in each side in !dim
	float get_hwire_spacing  () const {return 0.75*get_bar_extend();}
	float get_vwire_spacing  () const {return 0.25*get_bar_extend();}
	float get_standoff_height() const {return 0.3*pole_radius;}
	bool has_dim_set(unsigned d) const {return (dims & (1<<d));}
	cube_t calc_cbar(bool d) const;
public:
	power_pole_t(point const &base_, point const &center_, float pole_radius_, float height, float wires_offset_,
		float const pole_spacing_[2], uint8_t dims_, bool at_grid_edge_, bool const at_line_end_[2], bool residential_);
	bool is_at_grid_edge() const {return at_grid_edge;}
	bool has_transformer() const {return (dims == 3);}
	point get_top() const {return point(base.x, base.y, bcube.z2());}
	void get_wires_conn_pts(point pts[3], bool d) const;
	void get_top_wires_conn_pts(point pts[3]) const {get_wires_conn_pts(pts, ((dims&1) ? 0 : 1));} // use X if enabled, otherwise Y
	float get_pole_radius() const {return pole_radius;}
	float get_bsphere_radius(bool shadow_only) const {return (shadow_only ? radius : bsphere_radius);} // non-shadow pass includes wires bsphere radius
	cube_t const &get_outer_bcube() const {return bcube_with_wires;}
	cube_t get_bird_bcube  () const {return get_ped_occluder();} // centered on the pole; same as pedestrians
	cube_t get_ped_occluder() const;
	point get_nearest_connection_point(point const &to_pos, bool near_power_pole) const;
	void get_top_wire_end_pts(point top_wires[2][3][2]) const;
	point get_transformer_center() const;
	bool add_wire(point const &p1, point const &p2, bool add_pole);
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct hcap_space_t : public oriented_city_obj_t { // handicap space
	hcap_space_t(point const &pos_, float radius_, bool dim_, bool dir_);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};
struct hcap_with_dist_t : public hcap_space_t {
	float dmin_sq;
	hcap_with_dist_t(hcap_space_t const &hs, cube_t const &plot, vect_cube_t &bcubes, unsigned bcubes_end);
	bool operator<(hcap_with_dist_t const &v) const {return (dmin_sq < v.dmin_sq);}
};

struct manhole_t : public city_obj_t {
	manhole_t(point const &pos_, float radius_);
	float get_height() const {return 0.01*radius;}
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct mailbox_t : public model_city_obj_t {
	mailbox_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id()) {}
	float get_height() const {return 2.0*radius;}
	virtual unsigned get_model_id() const {return OBJ_MODEL_MAILBOX;}
};

struct bicycle_t : public model_city_obj_t {
	bicycle_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id()) {}
	virtual unsigned get_model_id() const {return OBJ_MODEL_BICYCLE;}
};

struct swingset_t : public model_city_obj_t {
	float anim_time=0.0, anim_scale=0.0;

	swingset_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id()) {}
	virtual unsigned get_model_id() const {return OBJ_MODEL_SWINGSET;}
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	void next_frame(point const &camera_bs);
};

struct trampoline_t : public model_city_obj_t {
	trampoline_t(point const &pos_, float height, rand_gen_t &rgen);
	virtual unsigned get_model_id() const {return OBJ_MODEL_TRAMPOLINE;}
};

struct dumpster_t : public model_city_obj_t {
	dumpster_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id()) {}
	virtual unsigned get_model_id() const {return OBJ_MODEL_DUMPSTER;}
};

struct umbrella_t : public model_city_obj_t { // large beach umbrella
	umbrella_t(point const &pos_, float height, bool dim_, bool dir_) : model_city_obj_t(pos_, height, dim_, dir_, get_model_id(), 1) {} // is_cylinder=1
	virtual unsigned get_model_id   () const {return OBJ_MODEL_BIG_UMBRELLA;}
	virtual float get_xy_coll_radius() const {return 0.33*model_city_obj_t::get_xy_coll_radius();} // smaller radius because we only collide with the vertical pole
};

struct potted_plant_t : public multi_model_city_obj_t {
	potted_plant_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_select) :
		multi_model_city_obj_t(pos_, height, dim_, dir_, OBJ_MODEL_PLANT, model_select, 1) {} // is_cylinder=1
};

struct flower_t : public multi_model_city_obj_t {
	flower_t(point const &pos_, float height, bool dim_, bool dir_, unsigned model_select) :
		multi_model_city_obj_t(pos_, height, dim_, dir_, OBJ_MODEL_FLOWER, model_select, 1) {min_alpha = 0.5;} // is_cylinder=1
};

struct traffic_cone_t : public city_obj_t {
	traffic_cone_t(point const &pos_, float radius_);
	float get_height() const {return 2.0*radius;}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct pond_t : public city_obj_t {
	pond_t(point const &pos_, float x_radius, float y_radius, float depth);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
	bool point_contains_xy(point const &p) const;
};

struct walkway_t : public oriented_city_obj_t, public walkway_base_t {
	colorRGBA map_mode_color;

	walkway_t(bldg_walkway_t const &w);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct pillar_t : public city_obj_t { // for walkway support
	bool is_concrete;

	pillar_t(cube_t const &bcube_, bool is_concrete_) : city_obj_t(bcube_), is_concrete(is_concrete_) {}
	float get_cylin_radius() const {return 0.25*(bcube.dx() + bcube.dy());}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct parking_solar_t : public oriented_city_obj_t {
	unsigned num_spaces, num_rows;
	parking_solar_t(cube_t const &c, bool dim_, bool dir_, unsigned ns, unsigned nr);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
	vect_cube_t const &get_legs() const;
};

struct city_bird_base_t : public city_obj_t {
	vector3d dir;
	city_bird_base_t(point const &pos_, float height, vector3d const &dir, unsigned model_id);
};

struct pigeon_t : public city_bird_base_t {
	pigeon_t(point const &pos_, float height, vector3d const &dir_) : city_bird_base_t(pos_, height, dir_, OBJ_MODEL_PIGEON) {}
	static void pre_draw (draw_state_t &dstate, bool shadow_only) {disable_hemi_lighting_pre_post(dstate, shadow_only, 0);}
	static void post_draw(draw_state_t &dstate, bool shadow_only) {disable_hemi_lighting_pre_post(dstate, shadow_only, 1);}
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

class city_obj_placer_t;
float const BIRD_ZVAL_ADJ = -0.05;

class city_bird_t : public city_bird_base_t {
	bool hit_min_alt=0;
	uint8_t state=0;
	unsigned loc_ix=0;
	float height=0.0, anim_time=0.0, takeoff_time=0.0, start_end_zmax=0.0, next_poop_time=0.0;
	vector3d velocity, dest_dir;
	point dest, prev_frame_pos;
	colorRGBA color;

	bool is_close_to_player() const;
	bool is_anim_cycle_complete(float new_anim_time) const;
	bool in_landing_dist() const;
	unsigned get_model_anim_id() const {return state;}
	void set_takeoff_time(rand_gen_t &rgen);
	void adjust_new_dest_zval();
	bool check_for_mid_flight_coll(float dir_dp, city_obj_placer_t &placer, rand_gen_t &rgen);
public:
	city_bird_t(point const &pos_, float height_, vector3d const &init_dir, colorRGBA const &color_, unsigned loc_ix_, rand_gen_t &rgen);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	void next_frame(float timestep, float delta_dir, point const &camera_bs, bool &tile_changed, bool &bird_moved, city_obj_placer_t &placer, rand_gen_t &rgen);
	bool dest_valid() const {return (dest != zero_vector);}
};

struct sign_t : public oriented_city_obj_t {
	bool two_sided, emissive, small, scrolling, free_standing;
	int sign_id=-1; // >= 0 means this is one of multiple alternatives
	colorRGBA bkg_color, text_color;
	cube_t connector, frame_bcube, text_bcube;
	string text;
	vector<float> char_pos; // used for scrolling while drawing

	sign_t(cube_t const &bcube_, bool dim_, bool dir_, string const &text_, colorRGBA const &bc, colorRGBA const &tc,
		bool two_sided_=0, bool emissive_=0, bool small_=0, bool scrolling_=0, bool fs=0, cube_t const &conn=cube_t());
	bool is_hospital_sign() const {return (text == "Hospital");}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	void draw_text(draw_state_t &dstate, city_draw_qbds_t &qbds, string const &text_to_draw, float first_char_clip_val=0.0, float last_char_clip_val=0.0) const;
};

struct stopsign_t : public oriented_city_obj_t {
	unsigned num_way;

	stopsign_t(point const &pos_, float height, float width, bool dim_, bool dir_, unsigned num_way_);
	cube_t get_bird_bcube() const;
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct city_flag_t : public oriented_city_obj_t {
	cube_t flag_bcube;
	point pole_base;
	float pole_radius;
	int flag_id; // -1 is default
	bool draw_as_model;

	city_flag_t(cube_t const &flag_bcube_, bool dim_, bool dir_, point const &pole_base_, float pradius, int flag_id_=-1);
	bool is_horizontal() const {return (flag_bcube.dz() > flag_bcube.get_sz_dim(!dim));}
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only) {disable_hemi_lighting_pre_post(dstate, shadow_only, 1);}
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct newsrack_t : public oriented_city_obj_t {
	colorRGBA color;
	unsigned style=0;

	newsrack_t(point const &pos_, float height, float width, float depth, bool dim_, bool dir_, unsigned style_, colorRGBA const &color_);
	cube_t get_bird_bcube() const;
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct park_path_t : public city_obj_t {
	vector<point> pts;
	float hwidth;
	colorRGBA color;
	cube_t plot;

	park_path_t(float hwidth_, colorRGBA const &color_, cube_t const &plot_) : hwidth(hwidth_), color(color_), plot(plot_) {}
	void calc_bcube_bsphere();
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool check_cube_coll_xy(cube_t const &c) const;
	bool check_point_contains_xy(point const &p) const;
};

class tile_drawer_t;

struct moving_walkway_t : public cube_t {
	bool dim, dir, active=1;
	float speed; // sped is always positive; dir is movement dir
	mutable int last_update_frame=0;
	mutable float move_time=0.0;
	cube_t track, sides[2];

	moving_walkway_t(cube_t const &c, bool dim_, bool dir_, float speed_=0.0);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, tile_drawer_t &td, bool shadow_only, bool draw_track, bool draw_sides) const;
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, point const &xlate, vector3d *cnorm) const;
};

struct skyway_t : public city_obj_t {
	bool valid=0, dim=0; // but no dir
	cube_t bot, top;
	vector<skyway_conn_t> ww_conns; // connection points to building walkways; ix envodes 2*dim + dir
	vect_cube_t entrances, sides, steps;
	vector<moving_walkway_t> mwws;

	skyway_t() {}
	void init(cube_t const &c, bool dim_);
	// Note: no pre_draw() and post_draw() because there can be only one
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
	void get_building_signs(vector<sign_t> &signs) const;
};

struct bird_place_t {
	point pos;
	vector3d orient;
	bool on_ground=0, in_use=0, use_orient=0;

	bird_place_t(point const &pos_, rand_gen_t &rgen) : pos(pos_), on_ground(1), use_orient(0) {set_rand_orient(rgen);} // ground constructor
	bird_place_t(point const &pos_, bool dim, bool dir, bool use_orient_=1) : pos(pos_), on_ground(0), use_orient(use_orient_)
	{orient[dim] = (dir ? 1.0 : -1.0);} // object constructor
	void set_rand_orient(rand_gen_t &rgen) {orient = rgen.signed_rand_vector_spherical();}
};
struct vect_bird_place_t : public vector<bird_place_t> {
	void add_placement(cube_t const &obj, bool dim, bool dir, bool orient_dir, float spacing, rand_gen_t &rgen);
	void add_placement_centerline(cube_t const &obj, bool dim, bool dir, rand_gen_t &rgen);
	void add_placement_rand_dim_dir(cube_t const &obj, float spacing, rand_gen_t &rgen);
	void add_placement_top_center(cube_t const &obj, rand_gen_t &rgen);
};

class bird_poop_manager_t {
	struct poop_t : public sphere_t {
		vector3d vel;
		poop_t(point const &pos_, float radius_, vector3d const &init_vel) : sphere_t(pos_, radius_), vel(init_vel) {}
	};
	struct splat_t : public sphere_t {
		float time;
		vector3d dir;
		splat_t(point const &pos_, float radius_, vector3d const &dir_) : sphere_t(pos_, radius_), time(tfticks), dir(dir_) {}	
		void add_quad(quad_batch_draw &qbd) const {qbd.add_quad_dirs(pos, radius*dir, radius*cross_product(dir, plus_z), WHITE, plus_z);}
	};
	vector<poop_t > poops;
	vector<splat_t> splats;
	cube_t city_bounds;
	quad_batch_draw splat_qbd;
	rand_gen_t rgen;

	void remove_oldest_splat();
public:
	void init(cube_t const &city_bounds_) {city_bounds = city_bounds_;}
	void add(point const &pos, float radius, vector3d const &init_vel) {poops.emplace_back(pos, radius, init_vel);}
	void next_frame();
	void draw(shader_t &s, vector3d const &xlate);
};

class city_obj_groups_t : public vector<cube_with_ix_t> {
	map<uint64_t, vector<unsigned>> by_tile;
	cube_t bcube;
public:
	cube_t const &get_bcube() const {return bcube;}
	void clear();
	void insert_obj_ix(cube_t const &c, unsigned ix);
	template<typename T> void add_obj(T const &obj, vector<T> &objs);
	template<typename T> void create_groups (vector<T>       &objs, cube_t &all_objs_bcube);
	template<typename T> void rebuild       (vector<T>       &objs, cube_t &all_objs_bcube);
	template<typename T> void update_obj_pos(vector<T> const &objs, cube_t &all_objs_bcube);
};

class city_obj_placer_t : private city_draw_qbds_t {
public: // road network needs access to parking lots and driveways for drawing
	vector<parking_lot_t> parking_lots;
	vector<driveway_t> driveways; // for houses
private:
	vector<bench_t> benches;
	vector<tree_planter_t> planters;
	vector<trashcan_t> trashcans;
	vector<fire_hydrant_t> fhydrants;
	vector<substation_t> sstations;
	vector<fountain_t> fountains;
	vector<divider_t> dividers; // dividers for residential plots
	vector<swimming_pool_t> pools;
	vector<pool_ladder_t> pladders;
	vector<pool_deck_t> pdecks;
	vector<power_pole_t> ppoles;
	vector<hcap_space_t> hcaps; // handicap signs painted on parking lots
	vector<manhole_t> manholes;
	vector<mailbox_t> mboxes;
	vector<traffic_cone_t> tcones;
	vector<pigeon_t> pigeons;
	vector<city_bird_t> birds;
	vector<sign_t> signs;
	vector<stopsign_t> stopsigns;
	vector<city_flag_t> flags;
	vector<newsrack_t> newsracks;
	vector<park_path_t> ppaths;
	vector<swingset_t> swings;
	vector<trampoline_t> tramps;
	vector<umbrella_t> umbrellas;
	vector<bicycle_t> bikes;
	vector<dumpster_t> dumpsters;
	vector<potted_plant_t> plants;
	vector<flower_t> flowers;
	vector<pond_t> ponds;
	vector<walkway_t> walkways;
	vector<pillar_t> pillars;
	vector<parking_solar_t> p_solars;
	// index is last obj in group
	city_obj_groups_t bench_groups, planter_groups, trashcan_groups, fhydrant_groups, sstation_groups, fountain_groups, divider_groups, pool_groups, plad_groups,
		pdeck_groups, ppole_groups, hcap_groups, manhole_groups, mbox_groups, tcone_groups, pigeon_groups, bird_groups, sign_groups, stopsign_groups, flag_groups,
		nrack_groups, ppath_groups, swing_groups, tramp_groups, umbrella_groups, bike_groups, dumpster_groups, plant_groups, flower_groups, pond_groups, walkway_groups,
		pillar_groups, p_solar_groups;
	skyway_t skyway; // optional
	bird_poop_manager_t bird_poop_manager;
	vector<city_zone_t> sub_plots; // reused across calls
	cube_t all_objs_bcube;
	vect_bird_place_t bird_locs;
	rand_gen_t bird_rgen;
	unsigned num_spaces=0, filled_spaces=0, num_x_plots=0, num_y_plots=0;
	float plot_subdiv_sz=0.0;
	bool has_residential_plots=0;
	
	bool gen_parking_lots_for_plot(cube_t const &full_plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix, vect_cube_t &bcubes, vect_cube_t &colliders, rand_gen_t &rgen);
	void add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen);
	void place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &tree_pos, rand_gen_t &rgen, unsigned buildings_end);
	void place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> const &tree_pos,
		vect_cube_t const &pond_blockers, rand_gen_t &rgen, bool have_streetlights);
	void place_residential_plot_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<road_t> const &roads,
		vect_cube_t const &pool_blockers, unsigned driveways_start, unsigned city_ix, rand_gen_t &rgen);
	bool place_swimming_pool(road_plot_t const &plot, city_zone_t const &yard, cube_with_ix_t const &house, bool dim, bool dir, bool shrink_dim, unsigned prev_blockers_end,
		float divider_hwidth, float const translate_dist[2], vect_cube_t const &pool_blockers, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen);
	bool check_bird_walkway_clearance(cube_t const &bc) const;
	void place_birds(cube_t const &city_bcube, rand_gen_t &rgen);
	void add_house_driveways(road_plot_t const &plot, vect_cube_t &temp_cubes, rand_gen_t &rgen, unsigned plot_ix);
	void place_stopsigns_in_isec(road_isec_t &isec);
	void place_objects_in_isec(road_isec_t &isec, bool is_residential, vector<point> const &hospital_signs, rand_gen_t &rgen);
	void add_ssign_and_slight_plot_colliders(vector<road_plot_t> const &plots, vector<road_isec_t> const isecs[3], vector<vect_cube_t> &plot_colliders) const;
	void add_objs_on_buildings(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &hospital_signs);
	template<typename T> void draw_objects(vector<T> const &objs, city_obj_groups_t const &groups, draw_state_t &dstate,
		float dist_scale, bool shadow_only, bool has_immediate_draw=0, bool draw_qbd_as_quads=0, float specular=0.75, float shininess=50.0);
	bool connect_power_to_point(point const &at_pos, bool near_power_pole);
	void connect_power_to_buildings(vector<road_plot_t> const &plots);
	bool check_walkway_coll_xy(point const &pos, float radius) const;
public:
	bool has_plot_dividers() const {return !dividers.empty();}
	bool have_animations  () const {return !birds   .empty();} // only birds are animated
	bool has_residential  () const {return has_residential_plots;}
	vector<power_pole_t> const &get_power_poles() const {return ppoles;} // used for city connectivity
	void set_plot_subdiv_sz(float sz) {plot_subdiv_sz = sz;}
	void gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars, vector<road_t> const &roads,
		vector<road_isec_t> isecs[3], cube_t const &city_bcube, unsigned city_id, bool have_cars, bool is_residential, bool have_streetlights);
	bool add_skyway(cube_t const &city_bcube, vect_bldg_walkway_t const &walkway_cands, rand_gen_t rgen);
	void finalize_streetlights_and_power(streetlights_t &sl, vector<vect_cube_t> &plot_colliders);
	static bool subdivide_plot_for_residential(cube_t const &plot, vector<road_t> const &roads,
		float plot_subdiv_sz, unsigned parent_plot_ix, unsigned city_ix, vect_city_zone_t &sub_plots);
	void draw_detail_objects(draw_state_t &dstate, bool shadow_only);
	void draw_transparent_objects(draw_state_t &dstate, bool shadow_only);
	bool proc_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
	bool get_color_at_xy(point const &pos, colorRGBA &color, bool skip_in_road) const;
	bool get_color_at_xy_pre_road(point const &pos, colorRGBA &color) const;
	void get_occluders(pos_dir_up const &pdu, vector3d const &xlate, vect_cube_t &occluders) const;
	void get_plot_cuts(cube_t const &plot, vect_cube_t &cuts) const;
	bool cube_int_underground_obj(cube_t const &c) const;
	bool move_to_not_intersect_driveway(point &pos, float radius, bool dim) const;
	void next_frame();
	void play_sounds();
	// birds
	int check_path_segment_coll(point const &p1, point const &p2, float radius) const;
	bool choose_bird_dest(point const &pos, float radius, unsigned &loc_ix, point &dest_pos, vector3d &dest_dir);
	void add_bird_poop(point const &pos, float radius, vector3d const &init_vel) {bird_poop_manager.add(pos, radius, init_vel);}
}; // city_obj_placer_t

float get_power_pole_offset();
bool check_city_building_line_coll_bs_any(point const &p1, point const &p2);

