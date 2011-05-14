// 3D World - Player game state classes
// by Frank Gennari
// 7/3/06

#ifndef _PLAYER_STATE_H_
#define _PLAYER_STATE_H_


#include "3DWorld.h"
using std::string;

int const NUM_WEAPONS = 15;

unsigned const POWERUP_TIME = unsigned(40*TICKS_PER_SECOND);


struct bbox { // size = 20

	float x1, y1, x2, y2;
	int index;
};


struct team_info { // size = 20

	bbox bb; // add others?
};


struct od_data { // size = 12

	int id, type, val;
	float dist;

	od_data() {};
	od_data(int type0, int id0, float dist0, int val0=0) : id(id0), type(type0), dist(dist0), val(val0) {};
	bool operator<(od_data const &o) const {return (dist < o.dist);}
};


struct waypoint_t {

	point pos;
	bool user_placed, visited;
	vector<float> smiley_times;
	vector<unsigned> next_wpts;

	waypoint_t(point const &p, bool const up);
	void mark_visited_by_smiley(unsigned const smiley_id);
	float get_time_since_last_visited(unsigned const smiley_id) const;
	void clear();
};


class waypt_used_set {

	unsigned last_wp;
	int last_frame;
	map<unsigned, int> used;

public:
	waypt_used_set() : last_wp(0), last_frame(0) {}	
	void clear();
	void insert(unsigned wp);
	bool is_valid(unsigned wp);
};


class unreachable_pts {

	int try_counts;
	float try_dist_sq;
	vector<point> cant_get;

public:
	unreachable_pts() : try_counts(0), try_dist_sq(0.0) {}
	
	void clear() {
		reset_try();
		cant_get.clear();
	}
	void reset_try(float tdist_sq=0.0) {
		try_counts  = 0;
		try_dist_sq = tdist_sq;
	}
	bool cant_reach(point const &pos) const;
	bool proc_target(point const &pos, point const &target, point const &last_target, bool can_reach);
	void shift_by(vector3d const &vd);
};


struct destination_marker {

	int xpos, ypos, dmin_sq;
	float min_depth;
	bool valid;

	destination_marker() : xpos(0), ypos(0), dmin_sq(0), min_depth(0.0) {}
	bool add_candidate(int x1, int y1, int x2, int y2, float depth, float radius);
	void update_dmin(int x, int y) {if (valid) dmin_sq = (xpos - x)*(xpos - x) + (ypos - y)*(ypos - y);}
	point get_pos() const;
	void clear() {valid = 0; dmin_sq = 0; min_depth = 0.0;}
};


struct player_state { // size = big

	bool plasma_loaded;
	int timer, target, objective, weapon, wmode, powerup, powerup_time, kills, deaths, cb_hurt, killer;
	int fire_frame, was_hit, hitter, target_visible, kill_time, rot_counter, uw_time;
	int ctv, target_type, stopped_time, last_waypoint;
	unsigned tid, fall_counter, chunk_index;
	float shields, plasma_size, zvel, dpos, last_dz, last_zvel;
	point target_pos, objective_pos, cb_pos;
	vector3d hit_dir, velocity;
	string name;
	int p_weapons[NUM_WEAPONS], p_ammo[NUM_WEAPONS];
	unsigned char *tdata;
	vector<int> balls;
	waypt_used_set waypts_used;
	unreachable_pts unreachable;
	destination_marker dest_mark;

	player_state() : plasma_loaded(0), timer(0), target(-1), objective(-1), weapon(0), wmode(0), powerup(0), powerup_time(0),
		kills(0), deaths(0), cb_hurt(0), killer(NO_SOURCE), fire_frame(0), was_hit(0), hitter(-2), target_visible(0),
		kill_time(0), rot_counter(0), uw_time(0), ctv(0), target_type(0), stopped_time(0), last_waypoint(-1), tid(0), fall_counter(0),
		chunk_index(0), shields(0.0), plasma_size(0.0), zvel(0.0), dpos(0.0), last_dz(0.0), last_zvel(0.0), target_pos(all_zeros),
		objective_pos(all_zeros), cb_pos(all_zeros), hit_dir(all_zeros), velocity(all_zeros), tdata(NULL) {}

	void init(bool w_start);
	bool no_weap() const;
	bool no_ammo() const;
	bool no_weap_or_ammo()   const {return (no_weap() || no_ammo());}
	float get_damage_scale() const {return ((powerup == PU_DAMAGE) ? 4.0 : 1.0);}
	float get_rspeed_scale() const {return ((powerup == PU_SPEED)  ? 1.5 : 1.0);}
	float get_fspeed_scale() const {return ((powerup == PU_SPEED)  ? 2.0 : 1.0);}
	float get_shield_scale() const {return ((powerup == PU_SHIELD) ? 0.5 : 1.0);}
};


#endif // _PLAYER_STATE_H_


