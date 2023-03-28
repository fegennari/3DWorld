// 3D World - physics ojbects classes
// by Frank Gennari
// 7/23/06
#pragma once

#include "3DWorld.h"
#include "collision_detect.h"

float const MAX_PART_CLOUD_RAD = 0.25;
float const DECAL_OFFSET       = 0.001;

extern float CAMERA_RADIUS, C_STEP_HEIGHT;

struct quad_batch_draw;


struct spark_t {
	float s=0.0;
	point pos;
	colorRGBA c=BLACK;
	static int  const status = 1;
	static float const radius;
	point const &get_pos() const {return pos;}

	spark_t(point const &p_, colorRGBA const &c_, float s_) : s(s_), pos(p_), c(c_) {}
	void draw(quad_batch_draw &qbd) const;
};


struct star { // size = 32
	float intensity=1.0;
	colorRGBA color;
	point pos;
};


struct obj_type { // size = 88
	int lifetime=0, tid=-1;
	unsigned flags=0;
	float mass=0, radius=0, volume=0, surface_area=0, density=0, air_factor=0, terminal_vel=0, friction_factor=0;
	float health=0, damage=0, min_t=0, max_t=0, elasticity=0, gravity=0, deform=0, def_recover=0;
	colorRGBA color;
};


struct basic_physics_obj { // size = 20
	point pos;
	int time;
	char status;

	basic_physics_obj(point const &pos_=all_zeros, int status_=0) : pos(pos_), time(0), status(status_) {}
	void init(point const &p);
	void init_gen_rand(point const &p, float rxy, float rz);
	point const &get_pos() const {return pos;}
	int get_replace_age() const {return time;}
	bool enabled() const {return (status != 0);}

	bool operator<(basic_physics_obj const &o) const {
		if (!enabled() &&  o.enabled()) return 0;
		if ( enabled() && !o.enabled()) return 1;
		return (o.time < time);
	}
};


struct bubble : public basic_physics_obj { // size = 44
	float radius=0.0, velocity=0.0;
	colorRGBA color;

	void gen(point const &p, float r=0.0, colorRGBA const &c=WATER_C);
	void draw(bool set_liquid_color) const;
	void apply_physics(unsigned i);
};

typedef vector<pair<float, unsigned> > order_vect_t;


struct particle_cloud : public basic_physics_obj { // size = 88
	struct part : public sphere_t {
		bool status=0;
	};
	bool acc_smoke=0, no_lighting=0, red_only=0;
	int source=-1, damage_type=0;
	float radius=0.0, init_radius=0.0, density=0.0, darkness=0.0, damage=0.0, timestep_factor=1.0;
	vector3d init_vel;
	colorRGBA base_color=BLACK;
	vector<part> parts;
	mutable vector<part> render_parts;
	static order_vect_t order;

	void gen(point const &p, colorRGBA const &bc, vector3d const &iv, float r, float den, float dark, float dam,
		int src, int dt, bool as, bool use_parts=1, bool nl=0, float spread=1.0, float tsfact=1.0);
	void draw(quad_batch_draw &qbd) const;
	void draw_part(point const &p, float r, colorRGBA c, quad_batch_draw &qbd) const;
	void apply_physics(unsigned i);
	float get_rscale() const {return CLIP_TO_01(1.0f - (radius - init_radius)/(MAX_PART_CLOUD_RAD - init_radius));}
	bool is_fire()     const {return (damage_type == BURNED || damage_type == FIRE);}
};


struct fire : public basic_physics_obj { // size = 60
	int source=-1;
	bool is_static=0;
	float radius=0.0, heat=0.0, cval=0.0, inten=0.0, light_bwidth=0.0;
	vector3d velocity;

	void gen(point const &p, float size, float intensity, int src, bool is_static_, float light_bw);
	void draw(quad_batch_draw &qbd, int &last_in_smoke) const;
	void apply_physics(unsigned i);
	void extinguish();
	int get_replace_age() const {return (is_static ? 0 : time);} // static fires are never replaced
};


class lightning_t { // size = 40
	struct lseg_t : public line3d {
		unsigned parent_len=0;
		float damage=0.0;
		bool full_path=1, hit_water=0, has_child=0;

		lseg_t() {color = LITN_C;}
		unsigned get_len() const {return (points.size() + parent_len);}
		void set_params(unsigned pl, float d, bool fp, bool hw, bool hc) {parent_len = pl; damage = d; full_path = fp; hit_water = hw; has_child = hc;}
	};
	bool enabled=0;
	int time=0;
	point hit_pos;
	vector<lseg_t> paths;
	rand_gen_t rgen;

	void gen_recur(point const &start, vector3d const &start_dir, float strength, unsigned parent_len);
public:
	void reset_time() {time = LITNING_TIME;}
	void disable();
	bool is_enabled() const {return (enabled && !paths.empty());}
	point const &get_hit_pos() const {return hit_pos;}
	void gen();
	void draw() const;
};


class physics_particle_manager {
protected:
	struct part_t { // size = 28
		point p; // position
		vector3d v; // velocity
		color_wrapper c;

		part_t() {}
		part_t(point const &p_, vector3d const &v_, colorRGBA const &c_) : p(p_), v(v_) {c.set_c4(c_);}
	};
	vector<part_t> parts;
public:
	void clear() {parts.clear();}
	void gen_particles(point const &pos, vector3d const &vadd, float vmag, float gen_radius, colorRGBA const &color, unsigned num);
	void apply_physics(float gravity, float terminal_velocity, bool emissive=0);
	void draw(float radius, int tid, bool emissive=0) const;
	bool is_pos_valid(point const &pos) const;
};


class water_particle_manager : public physics_particle_manager {
public:
	static colorRGBA calc_color(float mud_mix, float blood_mix) {
		return blend_color(BLOOD_C, blend_color(MUD_C, WATER_C, mud_mix, 1), blood_mix, 1);
	}
	void apply_physics();
	void draw() const;
};


struct decal_obj : public basic_physics_obj { // size = 76
	bool is_glass=0;
	int cid=-1, tid=-1, lifetime=0;
	float radius=0.0, alpha=1.0, rot_angle=0.0;
	colorRGBA color=BLACK;
	point ipos, cobj_cent_mass;
	vector3d orient;
	tex_range_t tex_range;

	void gen(point const &p, float r, float ang, vector3d const &o, int lt, int tid_, int cid_=-1, colorRGBA const &color_=BLACK,
		bool is_glass_=0, tex_range_t const &tr=tex_range_t());
	bool draw(quad_batch_draw &qbd) const;
	void maybe_draw_blood_trail(line_tquad_draw_t &blood_tqd) const;
	bool is_on_cobj(int cobj, vector3d *delta=nullptr) const;
	void check_cobj();
	void apply_physics(unsigned i);
	float get_alpha() const;
	vector3d get_platform_delta() const;
	point const get_pos() const {return (pos + get_platform_delta());}
};


struct dwobject : public basic_physics_obj { // size = 67(68) (dynamic world object)

	int coll_id=-1;
	short type=0, source=NO_SOURCE, flags=0;
	unsigned char direction=0;
	float health=0.0, angle=0.0;
	vector3d velocity, orientation, init_dir, vdeform;

	dwobject() : orientation(plus_z), init_dir(plus_z) {}
	dwobject(int type_, point const &pos_, vector3d const &vel_=all_zeros, int status_=0, float health_=0.0)
		: basic_physics_obj(pos_, status_), type(type_), health(health_), velocity(vel_), orientation(-plus_z), init_dir(-plus_z) {}
	float get_true_radius() const;
	float get_true_density() const;
	float get_true_mass() const;
	void advance_object(bool disable_motionless_objects, int iter, int obj_index);
	int surface_advance();
	void set_orient_for_coll(vector3d const *const forced_norm);
	int check_water_collision(float vz_old);
	void surf_collide_obj() const;
	void elastic_collision(point const &obj_pos, float energy, int obj_type);
	int object_bounce(int coll_type, vector3d &norm, float elasticity2, float z_offset, vector3d const &obj_vel=zero_vector);
	int object_still_stopped(int obj_index);
	void do_coll_damage();
	int check_vert_collision(int obj_index, int do_coll_funcs, int iter, vector3d *cnorm=NULL,
		vector3d const &mdir=all_zeros, bool skip_dynamic=0, bool only_drawn=0, int only_cobj=-1, bool skip_movable=0);
	int multistep_coll(point const &last_pos, int obj_index, unsigned nsteps);
	void update_vel_from_damage(vector3d const &dv);
	void damage_object(float damage, point const &dpos, point const &shoot_pos, int weapon);
	void update_precip_type();
	void add_obj_dynamic_light(int index) const;
	bool disabled() const {return (status == OBJ_STAT_BAD || status == OBJ_STAT_RES);}
	void disable() {status = 0; time = 0; health = 0.0;} // set all three just in case
	void verify_data() const {if (is_nan(pos) || is_nan(velocity)) print_and_terminate();}
	void print_and_terminate() const;
	bool lm_coll_invalid() const;
	bool invalid_coll(coll_obj const &cobj) const {return (type == LANDMINE && lm_coll_invalid() && cobj.is_player());}
	bool proc_stuck(bool static_top_coll);
	bool is_flat() const;
};


class vert_coll_detector {

	dwobject &obj;
	int type=0, iter=0;
	bool player=0, already_bounced=0, skip_dynamic=0, only_drawn=0, skip_movable=0;
	int coll=0, obj_index=0, do_coll_funcs=0, only_cobj=0;
	unsigned cdir=0, lcoll=0;
	float z_old=0.0, o_radius=0.0, z1=0.0, z2=0.0;
	point pos, pold;
	vector3d motion_dir, obj_vel;
	vector3d *cnorm;
	dwobject temp;

	bool safe_norm_div(float rad, float radius, vector3d &norm);
	void check_cobj_intersect(int index, bool enable_cfs, bool player_step);
	void init_reset_pos();
public:
	vert_coll_detector(dwobject &obj_, int obj_index_, int do_coll_funcs_, int iter_, vector3d *cnorm_,
		vector3d const &mdir=zero_vector, bool skip_dynamic_=0, bool only_drawn_=0, int only_cobj_=-1, bool skip_movable_=0) :
	obj(obj_), type(obj.type), iter(iter_), player(type == CAMERA || type == SMILEY || type == WAYPOINT), skip_dynamic(skip_dynamic_), only_drawn(only_drawn_),
		skip_movable(skip_movable_), obj_index(obj_index_), do_coll_funcs(do_coll_funcs_), only_cobj(only_cobj_), z_old(obj.pos.z), 
		pos(obj.pos), pold(obj.pos), motion_dir(mdir), obj_vel(obj.velocity), cnorm(cnorm_) {}

	void check_cobj(int index);
	int check_coll();
};


struct enabled_pos {
	point pos;
	bool enabled;

	enabled_pos() : enabled(0) {}
	enabled_pos(point const &pos_, bool enabled_=1) : pos(pos_), enabled(enabled_) {}
};


template<typename T> class obj_vector_t : public vector<T> {

	unsigned cur_avail=0;
	bool enabled=0;

	struct int_uint_pair {
		int i;
		unsigned u;
		int_uint_pair(int i_=0, unsigned u_=0) : i(i_), u(u_) {}
		bool operator<(int_uint_pair const &p) const {return (i > p.i);} // greater than on the int
	};
	void inc_cur_avail() {
		++cur_avail;
		if (cur_avail == size()) cur_avail = 0;
	}
public:
	using vector<T>::size;
	using vector<T>::empty;
	obj_vector_t(unsigned sz=0) : vector<T>(sz) {}

	unsigned choose_element(bool peek=0) {
		assert(!empty());
		assert(cur_avail < size());
		vector<T> const &v(*this);
		unsigned const start(cur_avail), sz((unsigned)size());

		for (unsigned i = 0; i < sz; ++i) {
			unsigned ix(i + start);
			if (ix >= sz) ix -= sz;
			if (!v[ix].enabled()) {cur_avail = ix; break;}
			if (v[ix].get_replace_age() > v[cur_avail].get_replace_age()) {cur_avail = ix;} // replace oldest element
		}
		unsigned const chosen(cur_avail);
		assert(chosen < size());
		if (!peek) {inc_cur_avail();}
		enabled = 1;
		return chosen;
	}
	bool choose_element_not_newer_than(int min_time, unsigned &chosen) {
		chosen = choose_element(0);
		return (!(*this)[chosen].enabled() || (*this)[chosen].get_replace_age() >= min_time);
	}
	void choose_elements(vector<unsigned> &ixs, unsigned num) {
		assert(num <= size());
		assert(cur_avail < size());
		ixs.clear();
		if (num == 0) return;
		ixs.reserve(num);
		vector<T> const &v(*this);
		unsigned const start(cur_avail), sz((unsigned)size());

		for (unsigned i = 0; i < sz; ++i) {
			unsigned ix(i + start);
			if (ix >= sz) ix -= sz;
			if (v[ix].enabled()) continue; // used element
			ixs.push_back(ix);
			
			if (ixs.size() == num) {
				cur_avail = ix;
				inc_cur_avail();
				return; // we're done
			}
		}
		unsigned const num_rem(num - (unsigned)ixs.size());
		vector<int_uint_pair> time_ixs;
		time_ixs.reserve(sz);

		for (unsigned i = 0; i < sz; ++i) {
			if (v[i].enabled()) {time_ixs.push_back(int_uint_pair(v[i].get_replace_age(), i));}
		}
		assert(num_rem <= time_ixs.size());
		sort(time_ixs.begin(), time_ixs.end()); // sort largest to smallest by time
		for (unsigned i = 0; i < num_rem; ++i) {ixs.push_back(time_ixs[i].u);}
		enabled = 1;
	}

	bool any_active() {
		if (!enabled) return 0; // fast cached case

		for (unsigned i = 0; i < size(); ++i) {
 			if (vector<T>::operator[](i).status) return 1;
		}
		enabled = 0; // record inactive state, will be reset to 1 in choose_element[s]()
		return 0;
	}
};


class cloud_manager_t : public obj_vector_t<particle_cloud> {
	unsigned cloud_tid=0, fbo_id=0, txsize=0, tysize=0;
	float frustum_z=0.0, last_xy_scale=0.0;
	mutable cube_t bcube;

	void set_red_only(bool val) {for (iterator i = begin(); i != end(); ++i) i->red_only = val;}
public:
	~cloud_manager_t() {free_textures();}
	void create_clouds();
	void update_lighting();
	cube_t get_bcube() const;
	float get_max_xy_extent() const;
	bool create_texture(bool force_recreate);
	void free_textures();
	void draw(bool no_update);
	bool is_inited()    const {return (txsize > 0);}
	float get_z_plane() const {return get_bcube().d[2][1];}
};


struct transform_data; // forward reference
typedef std::shared_ptr<transform_data> p_transform_data;


struct predef_obj { // size = 28
	point pos;
	int type=0, obj_used=-1;
	float regen_time=0.0, cur_time=0.0;

	predef_obj(point const &pos_=all_zeros, int type_=0, float rtime=0) : pos(pos_), type(type_), regen_time(rtime) {}
};


class obj_group { // size = 36

	obj_vector_t<dwobject> objects;
	vector<predef_obj> predef_objs;
	p_transform_data td;

public:
	unsigned init_objects=0, max_objs=0, app_rate=0, end_id=0, new_id=0;
	bool enabled=0, reorderable=0, predef_use_once=0;
	short type=0;
	unsigned char flags=0;

	void create(int obj_type_, unsigned max_objects_, unsigned init_objects_, unsigned app_rate_,
		bool init_enabled_, bool reorderable_, bool auto_max, bool predef_use_once_);
	unsigned get_updated_max_objs() const;
	void update_app_rate(float const val, unsigned min_app, unsigned max_app);
	void init_group();
	void preproc_this_frame();
	void reap_predef_objs(set<unsigned> const *const only_of_type=nullptr);
	void remove_reset_cobjs();
	size_t max_objects() const {return objects.size();}
	int choose_object(bool peek=0);
	void create_object_at(unsigned i, point const &pos);
	void enable();
	void disable();
	bool is_enabled()   const {return enabled;}
	bool large_radius() const;
	void set_enable(bool enabled_) {if (enabled_) enable(); else disable();}
	void toggle_enable() {set_enable(!enabled);}
	void free_objects();
	void shift(vector3d const &vd);
	p_transform_data get_td() {assert(td); return td;}
	dwobject       &get_obj(unsigned i)       {assert(enabled); assert(i < objects.size()); return objects[i];}
	dwobject const &get_obj(unsigned i) const {assert(enabled); assert(i < objects.size()); return objects[i];}
	bool obj_within_dist(unsigned i, point const &pos, float dist) const;
	bool temperature_ok() const;
	bool obj_has_shadow(unsigned obj_id) const;
	int get_ptype() const;
	void add_predef_obj(point const &pos, int type, int rtime);
	int get_next_predef_obj(dwobject &obj, unsigned ix);
	vector<predef_obj> const &get_predef_objs() const {return predef_objs;}
};


template<typename T> void reset_status(vector<T> &objs) {
	for (unsigned i = 0; i < objs.size(); ++i) {objs[i].status = 0;}
}


class reflective_cobjs_t {

	struct map_val_t {
		unsigned tid=0, tsize=0, last_size=0, faces_valid=0;
		cube_t bcube; // cached bcube to determine when the cobj is moved
	};
	map<unsigned, map_val_t> cobjs; // maps cid to tid
public:
	void add_cobj(unsigned cid);
	bool remove_cobj(unsigned cid);
	void clear() {free_textures(); cobjs.clear();}
	void mark_faces_invalid();
	void free_textures();
	void create_textures();
	unsigned get_tid_for_cid  (unsigned cid) const;
	unsigned get_tsize_for_cid(unsigned cid) const;
};


// group flag bits
#define WAS_ADVANCED     0x02
#define JUST_INIT        0x01
#define PRECIPITATION    0x04
#define APP_FROM_LT      0x08

// object type flag bits
#define SEMI_TRANSPARENT 0x01
#define BLEND            0x02
#define SPECULAR         0x04
#define LOW_SPECULAR     0x08
#define SELECTABLE       0x10
#define NO_FALL          0x20
#define FALL_EVERYWHERE  0x40
#define TAIL_WHEN_FALL   0x80
#define OBJ_IS_DROP      0x100
#define OBJ_ROLLS        0x200
#define DEFORMABLE       0x400
#define VERTEX_DEFORM    0x800
#define OBJ_EXPLODES     0x1000
#define EXPL_ON_COLL     0x2000
#define COLL_DESTROYS    0x4000
#define IS_PRECIP        0x8000
#define NO_WATER_DAMAGE  0x10000
#define OBJ_IS_FLAT      0x20000
#define OBJ_IS_CYLIN     0x40000
#define NO_COLL_DAMAGE   0x80000
#define OBJ_NON_SOLID    0x100000
#define OBJ_RAND_DIR_XY  0x200000

// object flag bits
#define XY_STOPPED       0x01
#define IS_CUBE_FLAG     0x02
#define Z_STOPPED        0x04
#define FLOATING         0x08
#define UNDERWATER       0x10
#define CAMERA_VIEW      0x20
#define IN_WATER         0x40
#define SHADOWED         0x80

// object flags (second byte)
#define TYPE_FLAG        0x0100
#define IS_ON_ICE        0x0200
#define STATIC_COBJ_COLL 0x0400
#define OBJ_COLLIDED     0x0800
#define PLATFORM_COLL    0x1000
#define USER_PLACED      0x2000
#define WAS_FIRED        0x2000 // Note: duplicate value of USER_PLACED flag
#define WAS_PUSHED       0x4000
#define FROZEN_FLAG      0x8000

unsigned const XYZ_STOPPED(XY_STOPPED | Z_STOPPED);
unsigned const ALL_COLL_STOPPED(XYZ_STOPPED | STATIC_COBJ_COLL);

