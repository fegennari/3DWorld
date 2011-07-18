// 3D World - collision detection/collision object classes
// by Frank Gennari
// 7/23/06

#ifndef _COLLISION_DETECT_H_
#define _COLLISION_DETECT_H_


#include "3DWorld.h"

typedef void (*collision_func)(int, int, vector3d const &, point const &, float, int);

// object/collision object types/status
enum {COLL_NULL      = 0, COLL_CUBE,     COLL_CYLINDER, COLL_SPHERE,  COLL_CYLINDER_ROT, COLL_POLYGON, COLL_INVALID};
enum {COLL_UNUSED    = 0, COLL_FREED,    COLL_PENDING,  COLL_STATIC,  COLL_DYNAMIC,      COLL_NEGATIVE};
enum {OBJ_STAT_BAD   = 0, OBJ_STAT_AIR,  OBJ_STAT_COLL, OBJ_STAT_GND, OBJ_STAT_STOP,     OBJ_STAT_RES};
enum {COBJ_LIT_FALSE = 0, COBJ_LIT_TRUE, COBJ_LIT_UNKNOWN};

unsigned const OBJ_CNT_REM_TJ = 1;

unsigned char const QD_TAG_QUAD     = 0x00;
unsigned char const QD_TAG_GLOBAL   = 0x01;
unsigned char const QD_TAG_TRIANGLE = 0x02;
unsigned char const QD_TAG_DLIST    = 0x03;
unsigned char const QD_TAG_TEXTURE  = 0x04;
unsigned char const QD_TAG_OLD      = 0x08;


struct quad_div {

	unsigned char dim, dir, tag, sbs;
	unsigned face;
	quad_div(unsigned dim_, unsigned dir_, unsigned char tag_, unsigned face_, unsigned sbs_) :
		dim(dim_), dir(dir_), tag(tag_), face(face_), sbs(sbs_)
	{
		assert(dim < 256 && dir < 256 && sbs < 256);
	}
	bool operator<(const quad_div &A) const {
		if (A.dim < dim ) return 0;
		if (A.dim > dim ) return 1;
		if (A.dir < dir ) return 0;
		if (A.dir > dir ) return 1;
		if (A.tag < tag ) return 0;
		if (A.tag > tag ) return 1;
		if (A.sbs < sbs ) return 0;
		if (A.sbs > sbs ) return 1;
		return (A.face < face);
	}
};


struct lv_val {

	unsigned status;
	unsigned short n0, n1;
	unsigned char *nvals;
	lv_val(unsigned status_=0, unsigned short n0_=0, unsigned short n1_=0, unsigned char *nvals_=NULL)
		: status(status_), n0(n0_), n1(n1_), nvals(nvals_) {}
};


unsigned char test_all_light_val(unsigned lighted, unsigned val);
bool check_lv_req(unsigned lighted, unsigned val1, unsigned val2, bool inv);



inline unsigned get_light_val(unsigned lighted, unsigned i) {

	return ((lighted & (0xF << (i<<2))) >> (i<<2));
}


inline void set_light_val(unsigned &lighted, unsigned val, unsigned i) {

	lighted &= ~(0xF << (i<<2)); // clear old bits
	lighted |=  (val << (i<<2)); // set new bits
}


inline bool is_partial_shadow(unsigned lighted) { // some lighted == 3

	return (test_all_light_val(lighted, 3) != 0);
}


inline bool require_either_lightval(unsigned lighted, unsigned val1, unsigned val2) {

	return check_lv_req(lighted, val1, val2, 0);
}


inline bool require_any_lightval(unsigned lighted, unsigned val1, unsigned val2) {

	return check_lv_req(lighted, val1, val2, 1);
}


//#define USE_HASHMAP // seems slower


#ifdef USE_HASHMAP
#include <hash_map>

struct quad_div_hash : public stdext::hash_compare<quad_div> {
	size_t operator()(quad_div const &qd) const {
		return (20011*size_t(qd.dim) + 9887*size_t(qd.dir) + 6121*size_t(qd.tag) + 14401*size_t(qd.face));
	}
	bool operator()(quad_div const &qd1, quad_div const &qd2) const {return (qd1 < qd2);}
};

typedef stdext::hash_map<quad_div, lv_val, quad_div_hash> lvmap; // doesn't seem to be significantly better
#else
typedef map<quad_div, lv_val> lvmap;
#endif


class obj_layer { // size = 60

public:
	bool draw, shadow, swap_txy;
	float elastic, tscale, specular, shine, tdx, tdy, refract_ix, light_atten;
	int tid;
	collision_func coll_func;
	colorRGBA color;

	obj_layer(float e=0.0, colorRGBA const &c=WHITE, bool d=0, const collision_func cf=NULL, int ti=-1,
		float ts=1.0, float spec=0.0, float shi=0.0, float tx=0.0, float ty=0.0, float ri=1.0, float la=0.0) :
		draw(d), shadow(1), swap_txy(0), elastic(e), tscale(ts), specular(spec), shine(shi),
		tdx(0.0), tdy(0.0), refract_ix(ri), light_atten(la), tid(ti), coll_func(cf), color(c) {}

	bool equal_params(const obj_layer &cobj) const {
		return (color == cobj.color && tid == cobj.tid && tscale == cobj.tscale && specular == cobj.specular &&
			shine == cobj.shine && draw == cobj.draw && elastic == cobj.elastic && tdx == cobj.tdx &&
			tdy == cobj.tdy && swap_txy == cobj.swap_txy && refract_ix == cobj.refract_ix && coll_func == cobj.coll_func);
	}
};


class cobj_params : public obj_layer { // size = 68

public:
	int cf_index;
	unsigned char surfs;
	bool is_dynamic, no_coll;
	//obj_layer *layer;

	cobj_params() : cf_index(-1), is_dynamic(0), no_coll(0), surfs(0) {}
	cobj_params(float e, colorRGBA const &c, bool d, bool id, const collision_func cf, int ci,
		int ti=-1, float ts=1.0, int s=0, float spec=0.0, float shi=0.0, bool nc=0, float tx=0.0, float ty=0.0) :
		obj_layer(e, c, d, cf, ti, ts, spec, shi, tx, ty), cf_index(ci), surfs(s), is_dynamic(id), no_coll(nc) {}
};


class coll_obj : public cube_t { // size = 248

public:
	cobj_params cp; // could store unique cps in a set of material properties to reduce memory requirements slightly
	char type, destroy, status, lighted;
	float radius, radius2, thickness, volume;
	int counter, id, platform_id, group_id, waypt_id;
	short npoints;
	unsigned char last_coll, coll_type;
	bool fixed, is_billboard;
	point points[N_COLL_POLY_PTS];
	vector3d norm;
	vector<vector<int> > sobjs;
	vector<int> occluders;
	set<int> shadow_depends;
	lvmap lightmap;

	coll_obj() : type(COLL_NULL), destroy(0), status(COLL_UNUSED), lighted(COBJ_LIT_UNKNOWN), radius(0.0), radius2(0.0), thickness(0.0),
		volume(0.0), counter(0), id(-1), platform_id(-1), group_id(-1), waypt_id(-1), npoints(0), last_coll(0), coll_type(0), fixed(0),
		is_billboard(0), norm(all_zeros) {}
	void init();
	void clear_lightmap(int mode, unsigned keep=0, bool keep_depends=0);
	void clear_lightmap_if_lighted_eq(int shadowed, int partial);
	void clear_dependent_cobjs_lightmaps(vector<coll_obj> &cobjs, unsigned ix) const;
	void update_shadowed_cobjs(vector<coll_obj> &cobjs, vector<int> const &indices, unsigned ix) const;
	void clear_internal_data(vector<coll_obj> &cobjs, vector<int> const &indices, unsigned ix);
	bool clear_lightmap_entry(lvmap::iterator it, int mode, unsigned keep, vector<pair<quad_div, lv_val> > *to_add=NULL);
	void calc_size();
	float calc_min_dim() const;
	bool clip_in_2d(float const bb[2][2], float &ztop, int d1, int d2, int dir) const;
	void set_npoints();
	void print_bounds() const;
	void bb_union(float bb[3][2], int init);
	void draw_cobj(unsigned i, int &last_tid, int &last_group_id, int &last_pri_dim);
	int  simple_draw(int ndiv, int in_cur_prim=PRIM_DISABLED, bool no_normals=0) const;
	void add_to_vector(vector<coll_obj> &cobjs, int type_);
	void check_if_cube();
	void add_as_fixed_cobj();
	int  add_coll_cobj();
	void re_add_coll_cobj(int index, int remove_old=1, int dhcm=0);
	void get_cvz_range(unsigned *zz, float zmin, float zmax, int x, int y) const;
	bool subdiv_fixed_cube(vector<coll_obj> &cobjs);
	int  intersects_cobj(coll_obj const &c, float toler=0.0) const;
	int  is_anchored() const;
	void shift_by(vector3d const &vd, bool force=0);
	void add_to_platform() const;
	bool dynamic_shadows_only() const;
	void add_shadow(unsigned light_sources, bool dynamic) const;
	bool cobj_plane_side_test(point const *pts, unsigned npts, point const &lpos) const;
	bool operator<(const coll_obj &cobj) const {return (volume < cobj.volume);} // sort by size
	bool equal_params(const coll_obj &c) const {return (type == c.type && status == c.status &&
		platform_id == c.platform_id && group_id == c.group_id && cp.equal_params(c.cp));}
	bool is_semi_trans()  const;
	bool no_draw()        const {return (!cp.draw || status == COLL_UNUSED || status == COLL_FREED || cp.color.alpha == 0.0);}
	bool disabled()       const {return (status != COLL_DYNAMIC && status != COLL_STATIC);}
	bool no_collision()   const {return (disabled() || cp.no_coll);}
	bool freed_unused()   const {return (status == COLL_FREED   || status == COLL_UNUSED);}
	bool all_shadowed()   const {return 0;} // *** WRITE - check lightmap for all 1s (cache result) ***
	bool is_occluder()    const {return (status == COLL_STATIC && type == COLL_CUBE && cp.draw && !is_semi_trans());}
	bool is_player()      const;
	bool is_invis_player()const;
	bool truly_static()   const;
	bool is_cylinder()    const {return (type == COLL_CYLINDER || type == COLL_CYLINDER_ROT);}
	bool can_be_scorched()const;
	point get_center_pt() const;
	float get_max_dim()   const;
	float get_light_transmit(point v1, point v2) const;
	void bounding_sphere(point &center, float &brad) const;
	bool has_poly_billboard_alpha() const;
	bool check_poly_billboard_alpha(point const &p1, point const &p2, float t) const;
	bool line_intersect(point const &p1, point const &p2) const;
	bool line_int_exact(point const &p1, point const &p2, float &t, vector3d &cnorm, float tmin, float tmax) const;
	void register_coll(unsigned char coll_time, unsigned char coll_type_) {last_coll = coll_time; coll_type = coll_type_;}

	// drawing code
	void draw_coll_cube(int do_fill, int tid) const;
	void draw_extruded_polygon(vector3d const *const normals, int tid) const;
	void draw_subdiv_cylinder(int nsides, int nstacks, bool draw_ends, bool no_bfc, bool no_lighting, int tid) const;
	void draw_subdiv_sphere_at(int ndiv, bool no_lighting, int tid) const;
};


unsigned const CLITE_FLAGS_DRAW = 0x0001; // draw
unsigned const CLITE_FLAGS_CUBE = 0x0002; // type cube
unsigned const CLITE_FLAGS_FIXD = 0x0004; // fixed
unsigned const CLITE_FLAGS_STAT = 0x0008; // static
unsigned const CLITE_FLAGS_DYNA = 0x0010; // dynamic
unsigned const CLITE_FLAGS_TRAN = 0x0020; // semi-transparent
unsigned const CLITE_FLAGS_CLER = 0x0040; // transparent (clear)
unsigned const CLITE_FLAGS_SMAL = 0x0080; // small volume
unsigned const CLITE_FLAGS_VSMA = 0x0100; // very small volume (< 0.0001)
unsigned const CLITE_FLAGS_DEST = 0x0200; // destroyable
unsigned const CLITE_FLAGS_VALD = 0x0400; // valid
unsigned const CLITE_FLAGS_NCOL = 0x0800; // no collisions


struct coll_cell { // size = 52

	float zmin, zmax, occ_zmin, occ_zmax;
	vector<int> cvals;
	vector<unsigned> cvz, indices;

	void clear(bool clear_vectors);
	void clear_cvz();
	void optimize(int x, int y);
	void update_opt(int x, int y) {if (!cvz.empty()) optimize(x, y);}
	void update_zmm(float zmin_, float zmax_, coll_obj const &cobj);
	
	void add_entry(int index) {
		if (INIT_CCELL_SIZE > 0 && cvals.capacity() == 0) cvals.reserve(INIT_CCELL_SIZE);
		cvals.push_back(index);
	}
};


class coll_cell_opt_batcher {

	bool enabled;
	set<pair<int, int> > to_proc;

public:
	void begin_batch();
	void end_batch();
	bool check_add_entry(int x, int y);
};


struct color_tid_vol : public cube_t {

	int cid, tid, destroy;
	bool draw, unanchored, is_2d;
	float volume, thickness, tscale;
	colorRGBA color;
	color_tid_vol(coll_obj const &cobj, float volume_, float thickness_, bool ua);
};


class platform { // animated (player controlled) scene object

	// constants
	bool const cont; // continuous - always in motion
	int shadow_mode;
	float const fspeed, rspeed; // velocity of forward/reverse motion in units per tick (can be negative)
	float const sdelay, rdelay; // start delay / reverse delay in ticks
	float const ext_dist, act_dist; // distance traveled, activation distance
	point origin; // initial position (moved by shifting)
	vector3d const dir; // direction of motion (travels to dir*dist)

	// state variables
	enum {ST_NOACT=0, ST_WAIT, ST_FWD, ST_CHDIR, ST_REV};
	int state; // 0 = not activated, 1 = activated but waiting, 2 = forward, 3 = waiting to go back, 4 = back
	float ns_time; // time to next state in ticks (can be negative if frame time is larger than a delay/travel time)
	point pos; // current position - dist is calculated from this point (delta = pos-origin)
	vector3d delta; // last change in position
	bool s_d_chg;

	// other data
	vector<unsigned> cobjs; // collision object(s) bound to this platform
	
public:
	platform(float fs=1.0, float rs=1.0, float sd=0.0, float rd=0.0, float dst=1.0, float ad=0.0,
		point const &o=all_zeros, vector3d const &dir_=plus_z, int sm=1, bool c=0);
	bool has_dynamic_shadows() const {return (shadow_mode >= 2 && (cont || state >= ST_FWD));}
	vector3d get_delta()       const {return (pos - origin);}
	vector3d get_range()       const {return dir*ext_dist;}
	vector3d get_last_delta()  const {return delta;}
	vector3d get_velocity()    const;
	void clear_last_delta() {delta = all_zeros;}
	void add_cobj(unsigned cobj);
	void next_frame();
	void shift_by(vector3d const &val);
	void reset();
	void activate();
	bool check_activate(point const &p, float radius);
	void advance_timestep();
	bool is_moving() const {return (state == ST_FWD || state == ST_REV);}
};


struct platform_cont : public deque<platform> {

	bool add_from_file(FILE *fp);
	void check_activate(point const &p, float radius);
	void shift_by(vector3d const &val);
	void advance_timestep();
};


class shadow_sphere : public sphere_t {

public:
	char lighted, ctype;
	int cid;

	shadow_sphere() : lighted(0), ctype(COLL_NULL), cid(-1) {}
	shadow_sphere(point const &pos0, float radius0, int cid0, bool lighted0);
	
	inline bool line_intersect(point const &p1, point const &p2) const {
		return (line_sphere_intersect(p1, p2, pos, radius) && (ctype == COLL_SPHERE || line_intersect_cobj(p1, p2)));
	}
	bool line_intersect_cobj(point const &p1, point const &p2) const;

	inline bool test_volume(point const *const pts, unsigned npts, point const &lpos) const {
		return (cid < 0 || test_volume_cobj(pts, npts, lpos));
	}
	bool test_volume_cobj(point const *const pts, unsigned npts, point const &lpos) const;
};


#endif // _COLLISION_DETECT_H_

