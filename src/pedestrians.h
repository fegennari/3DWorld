// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18
#pragma once

#include "3DWorld.h"

float const PED_WIDTH_SCALE  = 0.5; // ratio of collision radius to model radius (x/y)
float const PED_HEIGHT_SCALE = 2.5; // ratio of collision radius to model height (z)

extern double tfticks;

class ped_manager_t;

struct waiting_obj_t {
	float waiting_start;
	waiting_obj_t() : waiting_start(0.0) {}
	void reset_waiting() {waiting_start = tfticks;}
	float get_wait_time_secs() const {return (float(tfticks) - waiting_start)/TICKS_PER_SECOND;} // Note: only meaningful for cars stopped at lights or peds stopped at roads
};

// TODO: split into base class + city pedestrian vs. building person classes
struct pedestrian_t : public waiting_obj_t {

	point target_pos, dest_car_center; // since cars are sorted each frame, we can't find their positions by index so we need to cache them here
	vector3d dir, vel;
	point pos;
	float radius, speed, anim_time, retreat_time;
	unsigned plot, next_plot, dest_plot, dest_bldg; // Note: can probably be made unsigned short later, though these are global plot and building indices
	unsigned short city, model_id, ssn, colliding_ped, cur_rseed;
	unsigned char stuck_count;
	bool collided, ped_coll, is_stopped, in_the_road, at_crosswalk, at_dest, has_dest_bldg, has_dest_car, destroyed, in_building, following_player, is_on_stairs, has_key, is_female;

	pedestrian_t(float radius_) : target_pos(all_zeros), dir(zero_vector), vel(zero_vector), pos(all_zeros), radius(radius_), speed(0.0), anim_time(0.0), retreat_time(0.0),
		plot(0), next_plot(0), dest_plot(0), dest_bldg(0), city(0), model_id(0), ssn(0), colliding_ped(0), cur_rseed(1), stuck_count(0), collided(0), ped_coll(0), is_stopped(0),
		in_the_road(0), at_crosswalk(0), at_dest(0), has_dest_bldg(0), has_dest_car(0), destroyed(0), in_building(0), following_player(0), is_on_stairs(0), has_key(0), is_female(0) {}
	bool operator<(pedestrian_t const &ped) const {return ((city == ped.city) ? (plot < ped.plot) : (city < ped.city));} // currently only compares city + plot
	std::string get_name() const;
	std::string str() const;
	unsigned get_unique_id() const {return ssn;} // technically only unique if there are <= 65536 people
	float get_speed_mult() const;
	float get_height () const {return PED_HEIGHT_SCALE*radius;}
	float get_width  () const {return PED_WIDTH_SCALE *radius;}
	float get_coll_radius() const {return 0.6f*radius;} // using a smaller radius to allow peds to get close to each other
	float get_z1     () const {return (pos.z - radius);}
	float get_z2     () const {return (get_z1() + get_height());}
	cube_t get_bcube () const;
	point get_eye_pos() const {return (pos + vector3d(0.0, 0.0, (0.9*get_height() - radius)));}
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
	point get_dest_pos(cube_t const &plot_bcube, cube_t const &next_plot_bcube, ped_manager_t const &ped_mgr, int &debug_state) const;
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

