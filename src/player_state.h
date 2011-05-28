// 3D World - Player game state classes
// by Frank Gennari
// 7/3/06

#ifndef _PLAYER_STATE_H_
#define _PLAYER_STATE_H_


#include "3DWorld.h"
using std::string;

int const NUM_WEAPONS = 15;

unsigned const POWERUP_TIME = unsigned(40*TICKS_PER_SECOND);
float const S_SH_SCALE      = 2.0;


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


typedef unsigned short wpt_ix_t;
typedef vector<wpt_ix_t> waypt_adj_vect;


struct waypoint_t {

	bool user_placed, placed_item, goal, temp, visited;
	int came_from;
	float g_score, h_score, f_score;
	point pos;
	float last_smiley_time;
	waypt_adj_vect next_wpts, prev_wpts, visible_wpts;

	waypoint_t(point const &p=all_zeros, bool up=0, bool i=0, bool g=0, bool t=0);
	void mark_visited_by_smiley(unsigned const smiley_id);
	float get_time_since_last_visited(unsigned const smiley_id) const;
	void clear();
	bool unreachable() const {return prev_wpts.empty();}
	bool operator<(waypoint_t const &w) const {return (pos < w.pos);} // unused
};


struct wpt_goal {

	int mode; // 0: none, 1: user wpt, 2: placed item wpt, 3: goal wpt, 4: wpt index, 5: closest wpt, 6: closest visible wpt, 7: goal pos (new wpt)
	unsigned wpt;
	point pos;

	wpt_goal(int m=0, unsigned w=0, point const &p=all_zeros);
	bool is_reachable() const;
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
	void add(point const &pos) {cant_get.push_back(pos);}
	void shift_by(vector3d const &vd);
};


struct destination_marker {

	int xpos, ypos, dmin_sq;
	float min_depth;
	bool valid;

	destination_marker() : xpos(0), ypos(0), dmin_sq(0), min_depth(0.0), valid(0) {}
	bool add_candidate(int x1, int y1, int x2, int y2, float depth, float radius);
	void update_dmin(int x, int y) {if (valid) dmin_sq = (xpos - x)*(xpos - x) + (ypos - y)*(ypos - y);}
	point get_pos() const;
	void clear() {valid = 0; dmin_sq = 0; min_depth = 0.0;}
};


struct type_wt_t {
	unsigned type;
	float weight;
	type_wt_t(unsigned t=0, float w=1.0) : type(t), weight(w) {}
};


struct player_state { // size = big

	struct count_t {
		unsigned c;
		count_t(unsigned c_=0) : c(c_) {}
	};

	bool plasma_loaded, on_waypt_path;
	int timer, target, objective, weapon, wmode, powerup, powerup_time, kills, deaths, cb_hurt, killer;
	int init_frame, fire_frame, was_hit, hitter, target_visible, kill_time, rot_counter, uw_time;
	int target_type, stopped_time, last_waypoint;
	unsigned tid, fall_counter, chunk_index;
	float shields, plasma_size, zvel, dpos, last_dz, last_zvel, last_wpt_dist;
	point target_pos, objective_pos, cb_pos;
	vector3d hit_dir, velocity;
	string name;
	int p_weapons[NUM_WEAPONS], p_ammo[NUM_WEAPONS];
	unsigned char *tdata;
	vector<int> balls;
	map<unsigned, count_t> blocked_waypts;
	waypt_used_set waypts_used;
	unreachable_pts unreachable;
	destination_marker dest_mark;

	player_state() : plasma_loaded(0), on_waypt_path(0), timer(0), target(-1), objective(-1), weapon(0), wmode(0), powerup(0), powerup_time(0),
		kills(0), deaths(0), cb_hurt(0), killer(NO_SOURCE), init_frame(0), fire_frame(0), was_hit(0), hitter(-2), target_visible(0),
		kill_time(0), rot_counter(0), uw_time(0), target_type(0), stopped_time(0), last_waypoint(-1), tid(0), fall_counter(0),
		chunk_index(0), shields(0.0), plasma_size(0.0), zvel(0.0), dpos(0.0), last_dz(0.0), last_zvel(0.0), last_wpt_dist(0.0),
		target_pos(all_zeros), objective_pos(all_zeros), cb_pos(all_zeros), hit_dir(all_zeros), velocity(all_zeros), tdata(NULL) {}

	void init(bool w_start);
	void reset_wpt_state();
	bool no_weap() const;
	bool no_ammo() const;
	float weapon_range(bool use_far_clip) const;
	void verify_wmode();
	bool no_weap_or_ammo()   const {return (no_weap() || no_ammo());}
	float get_damage_scale() const {return ((powerup == PU_DAMAGE) ? 4.0 : 1.0);}
	float get_rspeed_scale() const {return ((powerup == PU_SPEED)  ? 1.5 : 1.0);}
	float get_fspeed_scale() const {return ((powerup == PU_SPEED)  ? 2.0 : 1.0);}
	float get_shield_scale() const {return ((powerup == PU_SHIELD) ? 0.5 : 1.0);}
	
	void smiley_fire_weapon(int smiley_id);
	int find_nearest_enemy(point const &pos, pos_dir_up const &pdu, point const &avoid_dir, int smiley_id,
		point &target, int &target_visible, float &min_dist) const;
	void check_cand_waypoint(point const &pos, point const &avoid_dir, int smiley_id,
		vector<od_data> &oddatav, unsigned i, int curw, float dmult, pos_dir_up const &pdu, bool next);
	int find_nearest_obj(point const &pos, pos_dir_up const &pdu, point const &avoid_dir, int smiley_id, point &target_pt,
		float &min_dist, vector<type_wt_t> types, int last_target_visible, int last_target_type);
	int check_smiley_status(dwobject &obj, int smiley_id);
	void drop_pack(point const &pos);
	int drop_weapon(vector3d const &coll_dir, vector3d const &nfront, point const &pos, int index, float energy, int type);
	void smiley_select_target(dwobject &obj, int smiley_id);
	float get_pos_cost(int smiley_id, point pos, point const &opos, pos_dir_up const &pdu, float radius, float step_height, bool check_dists);
	int smiley_motion(dwobject &obj, int smiley_id);
	void advance(dwobject &obj, int smiley_id);
	void shift(vector3d const &vd);
	void init_smiley_weapon(int smiley_id);
	bool target_in_range(point const &pos) const;
	void smiley_action(int smiley_id);

	// camera members
	void gamemode_fire_weapon();
	void switch_weapon(int val, int verbose);
	bool pickup_ball(int index);
	int fire_projectile(point fpos, vector3d dir, int shooter, int &chosen_obj);
	void update_camera_frame();
	void update_sstate_game_frame(int i);
	void free_balls();
	void update_weapon_cobjs(int i);
};


// function prototypes
bool check_step_dz(point &cur, point const &lpos, float radius);
int find_optimal_next_waypoint(unsigned cur, wpt_goal const &goal);
void find_optimal_waypoint(point const &pos, vector<od_data> &oddatav, wpt_goal const &goal);
bool can_make_progress(point const &pos, point const &opos);
bool is_valid_path(point const &start, point const &end);


#endif // _PLAYER_STATE_H_


