// 3D World - City Object Placement Header
// by Frank Gennari
// 08/14/21

#pragma once

#include "city.h"

struct city_draw_qbds_t {
	quad_batch_draw qbd, untex_qbd, untex_spec_qbd, emissive_qbd;
	bool empty() const {return (qbd.empty() && untex_qbd.empty() && untex_spec_qbd.empty() && emissive_qbd.empty());}
	bool has_untex_verts() const {return (!untex_qbd.empty() || !untex_spec_qbd.empty());}
};

struct city_obj_t : public sphere_t {
	cube_t bcube;
	city_obj_t() {}
	city_obj_t(point const &pos_, float radius_) : sphere_t(pos_, radius_) {}
	bool operator<(city_obj_t const &b) const {return (bcube.x1() < b.bcube.x1());} // sort by bcube x1
	cube_t const &get_outer_bcube() const {return bcube;}
	float get_bsphere_radius(bool shadow_only) const {return radius;}
	void set_bsphere_from_bcube() {*((sphere_t *)this) = bcube.get_bsphere();}
	static void post_draw(draw_state_t &dstate, bool shadow_only) {}
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const; // default, can be overridden in derived class
};

struct oriented_city_obj_t : public city_obj_t {
	bool dim, dir;
	oriented_city_obj_t(bool dim_=0, bool dir_=0) : dim(dim_), dir(dir_) {}
	oriented_city_obj_t(point const &pos_, float radius_, bool dim_=0, bool dir_=0) : city_obj_t(pos_, radius_), dim(dim_), dir(dir_) {}
};

struct bench_t : public oriented_city_obj_t {
	void calc_bcube();
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
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
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
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

struct substation_t : public oriented_city_obj_t {
	substation_t(cube_t const &bcube_, bool dim_, bool dir_);
	static void pre_draw (draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
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
	bool above_ground;

	swimming_pool_t(cube_t const &c, colorRGBA const &color_, colorRGBA const &wcolor_, bool above_ground_, bool dim_, bool dir_) :
		oriented_city_obj_t(c.get_cube_center(), c.get_bsphere_radius(), dim_, dir_), color(color_), wcolor(wcolor_), above_ground(above_ground_) {bcube = c;}
	float get_radius() const {assert(above_ground); return 0.25f*(bcube.dx() + bcube.dy());}
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
	point get_top() const {return point(base.x, base.y, bcube.z2());}
	void get_wires_conn_pts(point pts[3], bool d) const;
	void get_top_wires_conn_pts(point pts[3]) const {get_wires_conn_pts(pts, ((dims&1) ? 0 : 1));} // use X if enabled, otherwise Y
	float get_bsphere_radius(bool shadow_only) const {return (shadow_only ? radius : bsphere_radius);} // non-shadow pass includes wires bsphere radius
	cube_t const &get_outer_bcube() const {return bcube_with_wires;}
	cube_t get_ped_occluder() const;
	point get_nearest_connection_point(point const &to_pos, bool near_power_pole) const;
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

struct manhole_t : public city_obj_t {
	manhole_t(point const &pos_, float radius_);
	float get_height() const {return 0.01*radius;}
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct mailbox_t : public oriented_city_obj_t {
	mailbox_t(point const &pos_, float height, bool dim_, bool dir_);
	float get_height() const {return 2.0*radius;}
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct sign_t : public oriented_city_obj_t {
	bool two_sided, emissive;
	colorRGBA bkg_color, text_color;
	cube_t connector;
	string text;

	sign_t(cube_t const &bcube_, bool dim_, bool dir_, string const &text_, colorRGBA const &bc, colorRGBA const &tc, bool two_sided_=0, bool emissive_=0);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
};

struct city_flag_t : public oriented_city_obj_t {
	cube_t flag_bcube;
	point pole_base;
	float pole_radius;

	city_flag_t(cube_t const &flag_bcube_, bool dim_, bool dir_, point const &pole_base_, float pradius);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

class city_obj_groups_t : public vector<cube_with_ix_t> {
	map<uint64_t, vector<unsigned> > by_tile;
public:
	template<typename T> void add_obj(T const &obj, vector<T> &objs);
	template<typename T> void create_groups(vector<T> &objs, cube_t &all_objs_bcube);
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
	vector<divider_t> dividers; // dividers for residential plots
	vector<swimming_pool_t> pools;
	vector<power_pole_t> ppoles;
	vector<hcap_space_t> hcaps; // handicap signs painted on parking lots
	vector<manhole_t> manholes;
	vector<mailbox_t> mboxes;
	vector<sign_t> signs;
	vector<city_flag_t> flags;
	// index is last obj in group
	city_obj_groups_t bench_groups, planter_groups, trashcan_groups, fhydrant_groups, sstation_groups, divider_groups, pool_groups, ppole_groups,
		hcap_groups, manhole_groups, mbox_groups, sign_groups, flag_groups;
	vector<city_zone_t> sub_plots; // reused across calls
	cube_t all_objs_bcube;
	unsigned num_spaces, filled_spaces, num_x_plots, num_y_plots;
	float plot_subdiv_sz;
	
	struct cube_by_x1 {
		bool operator()(cube_t const &a, cube_t const &b) const {return (a.x1() < b.x1());}
	};
	bool gen_parking_lots_for_plot(cube_t plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix, vect_cube_t &bcubes, vect_cube_t &colliders, rand_gen_t &rgen);
	void add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen);
	void place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &tree_pos, rand_gen_t &rgen, unsigned buildings_end);
	void place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> const &tree_pos,
		rand_gen_t &rgen, bool is_residential, bool have_streetlights);
	void place_residential_plot_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, unsigned driveways_start, rand_gen_t &rgen);
	void add_house_driveways(road_plot_t const &plot, vect_cube_t &temp_cubes, rand_gen_t &rgen, unsigned plot_ix);
	void add_objs_on_buildings(unsigned city_id);
	template<typename T> void draw_objects(vector<T> const &objs, city_obj_groups_t const &groups, draw_state_t &dstate,
		float dist_scale, bool shadow_only, bool has_immediate_draw=0, bool draw_qbd_as_quads=0, float specular=0.75, float shininess=50.0);
	bool connect_power_to_point(point const &at_pos, bool near_power_pole);
	void connect_power_to_buildings(vector<road_plot_t> const &plots);
public:
	city_obj_placer_t() : num_spaces(0), filled_spaces(0), num_x_plots(0), num_y_plots(0), plot_subdiv_sz(0.0) {}
	bool has_plot_dividers() const {return !dividers.empty();}
	vector<power_pole_t> const &get_power_poles() const {return ppoles;} // used for city connectivity
	void clear();
	void set_plot_subdiv_sz(float sz) {plot_subdiv_sz = sz;}
	void gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars,
		unsigned city_id, bool have_cars, bool is_residential, bool have_streetlights);
	void move_and_connect_streetlights(streetlights_t &sl);
	static bool subdivide_plot_for_residential(cube_t const &plot, float plot_subdiv_sz, unsigned parent_plot_ix, vect_city_zone_t &sub_plots);
	void draw_detail_objects(draw_state_t &dstate, bool shadow_only);
	bool proc_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
	bool get_color_at_xy(point const &pos, colorRGBA &color, bool skip_in_road) const;
	void get_occluders(pos_dir_up const &pdu, vect_cube_t &occluders) const;
	void move_to_not_intersect_driveway(point &pos, float radius, bool dim) const;
};
