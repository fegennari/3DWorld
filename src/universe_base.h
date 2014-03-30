// 3D World - OpenGL CS184 Computer Graphics Project - base universe class definitions
// by Frank Gennari
// 9/4/02

#ifndef _UNIVERSE_BASE_H_
#define _UNIVERSE_BASE_H_


#include "3DWorld.h"


// special damage types
int const SWCLASS_UNDEF    = -1;
int const WCLASS_COLLISION = -2;
int const WCLASS_HEAT      = -3;
int const WCLASS_EXPLODE   = -4;
int const SCLASS_NONSHIP   = -5;
int const WCLASS_CAPTURE   = -6;
int const NO_OWNER         = -1;

float const USIZE_SCALE      = 4.0;
float const CELL_SIZE        = 100.0*USIZE_SCALE;
float const S_BODY_DENSITY   = 1.0E5;
float const EXPLOSION_DAMAGE = 50000.0;
float const MASS_SCALE       = 40000.0;
float const MAX_SOBJ_GRAVITY = 5.0;


// forward references
class free_obj;
class ship_coll_obj;


class cobj_vector_t : public vector<ship_coll_obj const *const> {

	void resize(size_t sz) {vector<ship_coll_obj const *const>::resize(sz, NULL);} // prohibited unless called from within this class
public:
	void clear();
	~cobj_vector_t() {clear();}
};


class uobject_base { // size = 28

public:
	float radius;
	upos_point_type pos; // required for high precision universe coordinates (can't use a sphere_t)

	uobject_base() : radius(0.0) {}
	upos_point_type const &get_pos() const {return pos;}
	float           get_radius()     const {return radius;}
	void set_pos(point const &pos_) {pos = pos_;}
};


class uobject : public uobject_base { // size = 44

public:
	char status;

	uobject() : status(0) {}
	virtual ~uobject() {}
	virtual void explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
		int align, unsigned eflags=0, free_obj const *parent_=NULL);
	void def_explode(float size, int etype, vector3d const &edir, int wclass=WCLASS_EXPLODE, int align=0,
		unsigned eflags=0, free_obj const *parent_=NULL, float dscale=1.0); // ALIGN_NEUTRAL==0 ?
	virtual float   get_bounding_radius() const {return radius;}
	bool            is_ok()               const {return (status != 1);}
	void add_gravity_vector_base(vector3d &vgravity, point const &mpos, float gfactor, float gmax) const;
	void gen_fragments(upos_point_type const &pos_offset=upos_point_type(0,0,0), float rscale=1.0) const;
	void gen_moving_fragments(point const &hit_pos, unsigned num, int tid=-1, float rscale=1.0, float vscale=1.0,
		vector3d const &vadd=zero_vector, colorRGBA const &pcolor=WHITE) const;
	float const *get_sphere_shadow_pmap(point const &sun_pos, point const &obj_pos, int ndiv) const;
	virtual std::string get_name() const = 0;
	virtual std::string get_info() const {return "";}
	virtual int get_id() const {return 0;} // standard uobject has no id
	virtual bool rename(std::string const &name_) {assert(0); return 0;}
	virtual int  get_type()  const {assert(0); return 0;}
	virtual int  get_owner() const {assert(0); return 0;}
	virtual cobj_vector_t const &get_cobjs() const;
	virtual bool casts_detailed_shadow() const {return !get_cobjs().empty();}
	virtual void draw_shadow_volumes(point const &targ_pos, float cur_radius, point const &sun_pos, int ndiv, bool test) const {assert(0);}
	virtual float get_radius_at(point const &p, bool exact=0) const {return radius;}
	virtual bool has_custom_shadow_profile() const {return 0;}
	virtual int get_fragment_tid(point const &hit_pos) const {return ROCK_SPHERE_TEX;}
	virtual bool sphere_intersection(point const &c, float r) const;
	virtual bool check_fragment_self_coll() const {return 0;}
};


struct ellipsoid_t {

	float xy_angle;
	vector3d scale, axis;
};


void reset_player_universe();
bool get_universe_sun_pos(point const &pos, point &spos);
bool has_sun_lighting(point const &pos);
int  set_uobj_color(point const &pos, float radius, bool known_shadowed, int shadow_thresh, point &sun_pos,
					uobject const *&sobj, float ambient_scale_s, float ambient_scale_no_s, shader_t *shader);
uobject *line_intersect_universe(point const &start, vector3d const &dir, float length, float line_radius, float &dist);


#endif


