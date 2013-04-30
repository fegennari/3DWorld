// 3D World - Ship models for universe mode - classes not in the free_obj class hierarchy
// by Frank Gennari
// 12/9/05

#ifndef _SHIP_UTIL_H_
#define _SHIP_UTIL_H_


#include "ship.h"
#include "obj_sort.h"


unsigned const BLOCK_SIZE = 1000;


extern unsigned alloced_fobjs[]; // testing


inline upos_point_type const &get_player_pos2() { // faster inlined version
	return player_ship().get_pos();
}

inline bool is_distant(point const &pos, float radius) {
	return (!dist_less_than(pos, player_ship().get_pos(), NDIV_SCALE_U*radius));
}

inline bool univ_sphere_vis_dist(point const &pos, float radius) {
	return (!is_distant(pos, radius) && univ_sphere_vis(pos, radius));
}



class usw_ray : public line_3dw {

	float w1, w2;
	colorRGBA color1, color2;

public:
	point prev, next;

	usw_ray() {}
	usw_ray(float w1_, float w2_, point const &p1_, point const &p2_, colorRGBA const &c1, colorRGBA const &c2)
		: line_3dw(p1_, p2_), w1(w1_), w2(w2_), color1(c1), color2(c2), prev(p1), next(p2)
	{
		assert(w1 > 0.0 && w2 > 0.0);
	}
	point const &get_pos() const {return p1;}
	void draw() const;
};


class us_fleet {

public:
	string name;
	unsigned align, ai, targ;
	float spread;
	point pos;
	u_ship *flagship;
	vector<pair<unsigned, unsigned> > ships; // sclass, count

	us_fleet() {}
	us_fleet(string const &name_, unsigned align_, unsigned ai_, unsigned targ_, float spread_,
		point const &pos_, unsigned counts[], unsigned multiplier=0);
	void set_flagship(unsigned sclass, float child_stray_dist);
	void spawn();
};


struct ship_cluster { // related to us_fleet, unused

	unsigned align;
	point center;
	vector3d velocity;
	vector<u_ship *> ships;

	void update();
};


struct query_data {
	
	vector<cached_obj> const *objs;
	point const pos;
	float const urm;
	float radius, damage, dist;
	unsigned eflags, index;
	int wclass, align;
	point ipos;
	free_obj const *parent, *fobj;
	uobject *ptr;
	bool exit_query;

	query_data(vector<cached_obj> const *const objs_, point const &pos_, float radius_, float urm_)
		: objs(objs_), pos(pos_), urm(urm_), radius(radius_), damage(0.0), dist(0.0),
		fobj(NULL), parent(NULL), ptr(NULL), exit_query(0) {}
};


struct base_query_data {

	vector<cached_obj> const *objs;
	point const pos; // compiler bug if const& ???
	free_obj const *const questioner;
	bool exit_query;

	base_query_data(vector<cached_obj> const *const objs_, point const &pos_, free_obj const *const questioner_) :
		objs(objs_), pos(pos_),	questioner(questioner_), exit_query(0) {}
};


struct closeness_data : public base_query_data {

	vector3d q_dir;
	float dmin, min_dist_sq, init_dmin, dscale;
	free_obj *closest;
	bool req_shields, req_dock, friendly;

	closeness_data(vector<cached_obj> const *const objs_, point const &pos_, float dmin_, float min_dist_sq_,
		free_obj const *const questioner_, bool req_sh=0, bool rdock=0, bool fr=0) :
		base_query_data(objs_, pos_, questioner_), q_dir(zero_vector), dmin(dmin_), init_dmin(dmin), min_dist_sq(min_dist_sq_),
		dscale(1.0), closest(NULL), req_shields(req_sh), req_dock(rdock), friendly(fr) {}
};


struct all_query_data : public base_query_data {

	float dmax_rscale, max_search_dist;
	vector<free_obj const*> &results;

	all_query_data(vector<cached_obj> const *const objs_, point const &pos_, float dmax, float search_dist,
		free_obj const *const questioner_, vector<free_obj const*> &results_) :
		base_query_data(objs_, pos_, questioner_), dmax_rscale(dmax), max_search_dist(dmax*search_dist), results(results_) {}
};


struct line_int_data {

	bool first_only, even_ncoll, visible_only, use_lpos;
	int check_parent;
	float length, line_radius, dist;
	point start, end, lpos;
	vector3d dir;
	uobject const *curr;
	free_obj const *ignore_obj;
	vector<uobject const*> *sobjs;

	line_int_data(point const &st, vector3d const &d, float len, uobject const *cur, free_obj const *ig, bool fo, int cp,
		float line_r=0.0) : first_only(fo), even_ncoll(0), visible_only(0), use_lpos(0), check_parent(cp), length(len),
		line_radius(line_r), dist(0.0), start(st), lpos(all_zeros), dir(d), curr(cur), ignore_obj(ig), sobjs(NULL) {}
};


template<typename T> class free_obj_block {

	T objs[BLOCK_SIZE];
	unsigned used, freed;
	bool in_use, valid;
public:
	free_obj_block() : used(0), freed(0), in_use(0), valid(1) {}

	T *alloc(unsigned type) {
		assert(valid);
		if (used == BLOCK_SIZE) return NULL;
		assert(used < BLOCK_SIZE);
		assert(type < T::max_type);
		assert(objs[used].status == 0);
		objs[used].set_type(type);
		if (VERIFY_REFS) objs[used].verify_status();
		return &objs[used++];
	}
	bool free_obj() {
		assert(valid);
		if (++freed == BLOCK_SIZE) {
			if (in_use) {
				for (unsigned i = 0; i < BLOCK_SIZE; ++i) objs[i].reset(); // is this necessary?
				used = freed = 0; // re-initialize
			}
			else {
				valid = 0; // in case someone tries to use this deleted block
				assert(alloced_fobjs[2] > 0);
				--alloced_fobjs[2];
				delete this; // is this legal?
			}
			return 0;
		}
		return 1;
	}
	void set_in_use(bool val) {assert(valid); in_use = val;}
};


template<typename T> class free_obj_allocator {

	free_obj_block<T> *last;
	free_obj_allocator(free_obj_allocator const &); // forbidden
	void operator=(free_obj_allocator const &); // forbidden
public:
	free_obj_allocator() : last(new free_obj_block<T>) {last->set_in_use(1);}

	T *alloc(unsigned type) {
		assert(last != NULL);
		T *ptr(last->alloc(type));
		
		if (ptr == NULL) {
			last->set_in_use(0); // unlock so that it can be freed
			last = new free_obj_block<T>;
			++alloced_fobjs[2];
			assert(last != NULL);
			last->set_in_use(1); // lock so that it can't be freed
			ptr = last->alloc(type);
			assert(ptr != NULL);
		}
		ptr->alloc_block = last;
		if (VERIFY_REFS) ptr->verify_status();
		return ptr;
	}

	~free_obj_allocator() {
	  //delete last; // what about the other pointers?
	}
};

// ship_config.cpp
void setup_ships();
u_ship *add_ship(unsigned sclass, unsigned align, unsigned ai, unsigned targ, point const &pos, float spread);

// universe_control.cpp
void send_warning_message(string const &msg);
void disable_player_ship();
void destroy_player_ship(bool captured);
bool rename_obj(uobject *obj, unsigned alignment);
uobject const *get_closest_world_ptr(point const &pos, int type);
uobject const *choose_dest_world(point const &pos, int exclude_id, unsigned align);
bool check_dest_ownership(int uobj_id, point const &pos, free_obj *own, bool check_for_land, bool homeworld);
void fire_planet_killer(u_ship const *const ship, point const &ship_pos, vector3d const &fire_dir, float fire_range, int obj_types);
orbiting_ship *add_orbiting_ship(unsigned sclass, bool guardian, bool on_surface, bool pos_from_parent, free_obj const *parent, urev_body *obj);

// ship.cpp
void print_n_spaces(int n);
u_ship *create_ship(unsigned sclass, point const &pos0, unsigned align, unsigned ai_type, unsigned target_mode, bool rand_orient);
bool add_uobj_ship(u_ship *ship);
bool add_uobj(free_obj *obj, int coll_test=0);
void reset_player_ship();
free_obj *line_intersect_free_objects(line_int_data &li_data, int obj_types, unsigned align, bool align_only);
uobject *line_intersect_objects(line_int_data &li_data, free_obj *&fobj, int obj_types);
unsigned check_for_obj_coll(point const &pos, float radius);
void get_all_close_objects(all_query_data &qdata);
void register_attack_from(free_obj const *attacker, unsigned target_align);
void register_damage(int t_sclass, int s_sclass, int wclass, float damage, unsigned s_align, unsigned t_align, bool is_kill);
void change_speed_mode(int val);
void toggle_player_ship_stop();
void change_fire_primary();
void reset_player_target();
void toggle_autopilot();
void toggle_dock_fighters();
void toggle_hold_fighters();
uparticle *gen_particle(unsigned type, colorRGBA const &c1, colorRGBA const &c2, unsigned lt, point const &pos,
						vector3d const &vel, float size, float damage, unsigned align, bool coll, int texture_id=-1);
void add_parts_projs(point const &pos, float radius, vector3d const &dir, colorRGBA const &color, int type, int align, free_obj const *const parent);
us_projectile *create_projectile(unsigned type, free_obj const *const parent, unsigned align, point const &pos,
								 vector3d const &vel, vector3d const &dir, vector3d const &upv);
void apply_explosion(point const &pos, float radius, float damage, unsigned eflags, int wclass, uobject *ptr, free_obj const *parent);
free_obj const *check_for_incoming_proj(point const &pos, int align, float dist);
void shift_univ_objs(point const &pos, bool shift_player_ship);
void draw_univ_objects(point const &pos);
void purge_old_objs();
void add_other_ships(int align, unsigned num=0, bool initial=0);
void gen_lightning_from(point const &pos, float radius, float dist, free_obj const *src);
void add_colored_lights(point const &pos, float radius, colorRGBA const &color, float time, unsigned num, free_obj const *const obj);

void add_br_light(unsigned index, point const &pos, float radius, free_obj const *const parent);
void apply_explosions();
void check_explosion_refs();

us_class const &get_sclass_obj(unsigned sclass);
float get_ship_cost     (unsigned sclass, unsigned align, unsigned reserve_credits=0, float discount=0.0);
bool alloc_resources_for(unsigned sclass, unsigned align, unsigned reserve_credits=0, float discount=0.0);

// draw_ship.cpp
void setup_colors_draw_flare(point const &pos, point const &xlate, float xsize, float ysize, colorRGBA const &color, int flare_tex=BLUR_TEX);
void draw_crosshair(upos_point_type const &pos, float dist, colorRGBA const &color);
void draw_crosshair_from_camera(point const &pos, colorRGBA const &color);
void add_lightning_wray(float width, point const &p1, point const &p2);


#endif

