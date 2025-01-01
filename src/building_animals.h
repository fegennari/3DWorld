// 3D World - Buildings Animal Class Definitions
// by Frank Gennari
// 3-22-23
#pragma once

#include "3DWorld.h"

template<typename T> cube_t get_cube_height_radius(point const &center, T radius, float height);


struct building_animal_t {
	point pos, last_pos;
	vector3d dir;
	float radius=0.0, speed=0.0, anim_time=0.0, wake_time=0.0, dist_since_sleep=0.0;
	unsigned id=0;
	bool shadow_non_visible=0;

	building_animal_t(float xval) : pos(xval, 0.0, 0.0) {}
	building_animal_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_) : pos(pos_), dir(dir_), radius(radius_), id(id_) {}
	bool operator<(building_animal_t const &a) const {return (pos.x < a.pos.x);} // compare only xvals
	bool is_moving  () const {return (speed     > 0.0);}
	bool is_sleeping() const {return (wake_time > 0.0);}
	vector3d get_upv() const {return plus_z;}
	void sleep_for(float time_secs_min, float time_secs_max);
	float move(float timestep, bool can_move_forward=1);
	bool detailed_sphere_coll(point const &sc, float sr, point &coll_pos, float &coll_radius) const {return 1;} // defaults to true
};

struct rat_t : public building_animal_t {
	point dest, fear_pos;
	float height=0.0, hwidth=0.0, fear=0.0;
	unsigned tunnel_tank_ix=0; // for sewer or pet store rats; could also add room_ix if it helps
	bool is_hiding=0, near_player=0, attacking=0, dead=0;

	// this first constructor is for the lower_bound() call in vect_rat_t::get_first_rat_with_x2_gt()
	rat_t(float xval) : building_animal_t(xval) {}
	rat_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_, bool dead_=0, unsigned tix=0);
	bool operator<(rat_t const &r) const {return (pos.x < r.pos.x);} // compare only xvals
	static bool allow_in_attic () {return 1;}
	static bool not_by_ext_door() {return 0;}
	float get_hlength() const {return radius;} // this is the bounding radius, so it represents the longest dim (half length)
	float get_height () const {return height;}
	float get_xy_radius() const {return radius;}
	point get_center () const {return point(pos.x, pos.y, (pos.z + 0.5f*height));}
	cube_t get_bcube () const {return get_cube_height_radius(pos, radius, height);} // used for collision detection and VFC; bounding cube across rotations
	cube_t get_bcube_with_dir() const; // used for model drawing; must be correct aspect ratio
	bool is_facing_dest() const;
};

struct spider_t : public building_animal_t {
	vector3d upv;
	point last_valid_pos;
	float update_time=0.0, web_start_zval=0.0, jump_vel_z=0.0, jump_dist=0.0;
	unsigned tunnel_tank_ix=0; // for sewer or pet store spiders
	bool on_web=0, web_dir=0, squished=0, in_tank=0; // web_dir: 0=going down, 1=going up

	// this first constructor is for the lower_bound() call in vect_rat_t::get_first_rat_with_x2_gt()
	spider_t(float xval) : building_animal_t(xval) {}
	spider_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_, unsigned tix=0, bool in_tank_=0);
	static bool allow_in_attic () {return 1;}
	static bool not_by_ext_door() {return 1;}
	float get_xy_radius() const {return 2.0*radius;}
	float get_height   () const {return 2.0*radius;}
	vector3d get_size  () const;
	cube_t get_bcube   () const; // used for collision detection and VFC
	vector3d get_upv   () const {return upv;}
	void choose_new_dir(rand_gen_t &rgen);
	// jumping logic
	void jump(float vel);
	bool is_jumping() const {return (jump_vel_z != 0.0);}
	void end_jump  ();
};

struct snake_t : public building_animal_t {
	// for snakes, pos is (xc, yc, z1), radius is the max body radius, and dir is the head direction and direction of movement
	float length=0.0, xy_radius=0.0;
	unsigned stuck_counter=0;
	bool has_rattle=0;
	vector3d last_valid_dir;
	colorRGBA color;
	vector<point> segments; // segment centers: first = head, last = tail

	snake_t(float xval) : building_animal_t(xval) {}
	snake_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_);
	static bool allow_in_attic () {return 0;}
	static bool not_by_ext_door() {return 1;}
	void  calc_xy_radius();
	float get_xy_radius () const {return xy_radius;} // must be fast
	float get_height    () const {return radius;}
	float get_seg_length() const {return length/segments.size();}
	float get_seg_radius(float seg_ix) const;
	point const &get_head_pos() const {assert(!segments.empty()); return segments.front();}
	point       &get_head_pos()       {assert(!segments.empty()); return segments.front();}
	cube_t get_bcube    () const;
	void move_segments(float dist);
	bool check_line_int_xy(point const &p1, point const &p2, bool skip_head, vector3d *seg_dir=nullptr) const;
	bool check_sphere_int    (point const &sc, float sr, bool skip_head, vector3d *seg_dir=nullptr, point *closest_pos=nullptr) const;
	bool detailed_sphere_coll(point const &sc, float sr, point &coll_pos, float &coll_radius) const;
	float get_curve_factor() const;
};

enum {INSECT_TYPE_FLY=0, INSECT_TYPE_ROACH, INSECT_TYPE_CENTIPEDE, NUM_INSECT_TYPES};

struct insect_t : public building_animal_t {
	vector3d delta_dir;
	float accel=0.0; // for flies
	float dist_to_sleep=0.0; // for roaches
	unsigned char type=0, stuck_counter=0;
	bool has_target=0, target_player=0; // for flies
	bool is_scared=0, no_scare=0, squished=0; // for roaches

	insect_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_, unsigned char type_=INSECT_TYPE_FLY) :
		building_animal_t(pos_, radius_, dir_, id_), type(type_) {}
	static bool allow_in_attic () {return 0;} // could allow it for flies but not cockroaches?
	static bool not_by_ext_door() {return 0;}
	bool flies() const {return (type == INSECT_TYPE_FLY);}
	float get_xy_radius() const {return radius;}
	float get_height   () const;
	float get_z2       () const {return (pos.z + 0.5*get_height());}
	vector3d get_orient() const {return vector3d(dir.x, dir.y, 0.0).get_norm();} // XY plane for all insects
	cube_t get_bcube () const; // used for collision detection and VFC; bounding cube across rotations
	cube_t get_bcube_with_dir() const; // used for model drawing; must be correct aspect ratio
};

template<typename T> struct vect_animal_t : public vector<T> {
	bool placed;
	float max_radius, max_xmove;
	vect_animal_t() : placed(0), max_radius(0.0), max_xmove(0.0) {}
	
	void add(T const &animal) {
		this->push_back(animal);
		this->back().id = this->size(); // rat_id starts at 1
		max_eq(max_radius, animal.get_xy_radius());
	}
	void do_sort() {
		sort(this->begin(), this->end()); // sort by xval
		max_xmove = 0.0; // reset for this frame
	}
	typename vector<T>::const_iterator get_first_with_xv_gt(float x) const {return std::lower_bound(this->begin(), this->end(), T(x));}
	void update_delta_sum_for_animal_coll(point const &pos, point const &cur_obj_pos, float radius,
		float z1, float z2, float radius_scale, float &max_overlap, vector3d &delta_sum) const;
};
typedef vect_animal_t<rat_t   > vect_rat_t;
typedef vect_animal_t<spider_t> vect_spider_t;
typedef vect_animal_t<snake_t > vect_snake_t;
typedef vect_animal_t<insect_t> vect_insect_t;

