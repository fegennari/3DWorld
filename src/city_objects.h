// 3D World - City Object Placement Header
// by Frank Gennari
// 08/14/21

#pragma once

#include "city.h"

struct city_obj_t : public sphere_t {
	cube_t bcube;
	city_obj_t() {}
	city_obj_t(point const &pos_, float radius_) : sphere_t(pos_, radius_) {}
	bool operator<(city_obj_t const &b) const {return (bcube.x1() < b.bcube.x1());} // sort by bcube x1
	static void post_draw(draw_state_t &dstate, bool shadow_only) {}
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const; // default, can be overridden in derived class
};

struct bench_t : public city_obj_t {
	bool dim, dir;
	bench_t() : dim(0), dir(0) {}
	void calc_bcube();
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const;
};

struct tree_planter_t : public city_obj_t {
	tree_planter_t(point const &pos_, float radius_, float height);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const;
};

struct fire_hydrant_t : public city_obj_t {
	float cylin_radius;
	vector3d orient;
	fire_hydrant_t(point const &pos_, float radius_, float height, vector3d const &orient_);
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	static void post_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const;
	bool proc_sphere_coll(point &pos_, point const &p_last, float radius_, point const &xlate, vector3d *cnorm) const;
};

struct divider_t : public city_obj_t {
	unsigned type;
	bool dim, dir;
	uint8_t skip_dims;
	divider_t(cube_t const &c, unsigned type_, bool dim_, bool dir_, unsigned sd=0) :
		city_obj_t(c.get_cube_center(), c.get_bsphere_radius()), type(type_), dim(dim_), dir(dir_), skip_dims(sd) {bcube = c;}
	static void pre_draw(draw_state_t &dstate, bool shadow_only);
	void draw(draw_state_t &dstate, quad_batch_draw &qbd, float dist_scale, bool shadow_only) const;
};

class city_obj_placer_t {
public: // road network needs access to parking lots and driveways for drawing
	vector<parking_lot_t> parking_lots;
	vector<driveway_t> driveways; // for houses
private:
	vector<bench_t> benches;
	vector<tree_planter_t> planters;
	vector<fire_hydrant_t> fire_hydrants;
	vector<divider_t> dividers; // dividers for residential plots
	vector<cube_with_ix_t> bench_groups, planter_groups, fire_hydrant_groups, divider_groups; // index is last object in group
	quad_batch_draw qbd;
	vect_cube_t temp_blockers;
	unsigned num_spaces, filled_spaces;
	
	struct cube_by_x1 {
		bool operator()(cube_t const &a, cube_t const &b) const {return (a.x1() < b.x1());}
	};
	bool gen_parking_lots_for_plot(cube_t plot, vector<car_t> &cars, unsigned city_id, unsigned plot_ix, vect_cube_t &bcubes, vect_cube_t &colliders, rand_gen_t &rgen);
	void add_cars_to_driveways(vector<car_t> &cars, vector<road_plot_t> const &plots, vector<vect_cube_t> &plot_colliders, unsigned city_id, rand_gen_t &rgen) const;
	void place_trees_in_plot(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> &tree_pos, rand_gen_t &rgen, unsigned buildings_end);
	template<typename T> void add_obj_to_group(T const &obj, cube_t const &bcube, vector<T> &objs, vector<cube_with_ix_t> &groups, bool &is_new_tile);
	template<typename T> void sort_grouped_objects(vector<T> &objs, vector<cube_with_ix_t> const &groups);
	void place_detail_objects(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, vector<point> const &tree_pos, rand_gen_t &rgen, bool is_new_tile, bool is_residential);
	void place_plot_dividers(road_plot_t const &plot, vect_cube_t &blockers, vect_cube_t &colliders, rand_gen_t &rgen, bool is_new_tile, float plot_subdiv_sz);
	void add_house_driveways(road_plot_t const &plot, vect_cube_t &temp_cubes, rand_gen_t &rgen, unsigned plot_ix);
	template<typename T> void draw_objects(vector<T> const &objs, vector<cube_with_ix_t> const &groups,
		draw_state_t &dstate, float dist_scale, bool shadow_only, bool not_using_qbd=0);
public:
	void clear();
	void gen_parking_and_place_objects(vector<road_plot_t> &plots, vector<vect_cube_t> &plot_colliders, vector<car_t> &cars,
		unsigned city_id, bool have_cars, bool is_residential, float plot_subdiv_sz);
	static bool subdivide_plot_for_residential(cube_t const &plot, float plot_subdiv_sz, vect_city_zone_t &sub_plots);
	void draw_detail_objects(draw_state_t &dstate, bool shadow_only);
	bool proc_sphere_coll(point &pos, point const &p_last, float radius, vector3d *cnorm) const;
	bool line_intersect(point const &p1, point const &p2, float &t) const;
	bool get_color_at_xy(point const &pos, colorRGBA &color) const;
};

inline uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}
