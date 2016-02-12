// 3D World - collision detection/collision object classes
// by Frank Gennari
// 7/23/06

#ifndef _COLLISION_DETECT_H_
#define _COLLISION_DETECT_H_


#include "function_registry.h"
#include "trigger.h"

typedef bool (*collision_func)(int, int, vector3d const &, point const &, float, int);

// object/collision object types/status
enum {COLL_NULL     = 0, COLL_CUBE,     COLL_CYLINDER, COLL_SPHERE,  COLL_CYLINDER_ROT, COLL_POLYGON,  COLL_CAPSULE, COLL_TORUS, COLL_INVALID};
enum {COLL_UNUSED   = 0, COLL_FREED,    COLL_PENDING,  COLL_STATIC,  COLL_DYNAMIC,      COLL_NEGATIVE, COLL_TO_REMOVE};
enum {OBJ_STAT_BAD  = 0, OBJ_STAT_AIR,  OBJ_STAT_COLL, OBJ_STAT_GND, OBJ_STAT_STOP,     OBJ_STAT_RES};
enum {COBJ_TYPE_STD = 0, COBJ_TYPE_MODEL3D, COBJ_TYPE_VOX_TERRAIN};

// collision object destroyability
enum {NON_DEST=0, DESTROYABLE, SHATTERABLE, SHATTER_TO_PORTAL, EXPLODEABLE};

unsigned const OBJ_CNT_REM_TJ = 1;

struct geom_xform_t;


struct base_mat_t { // size = 36
	
	int tid;
	float shine;
	colorRGBA color;
	colorRGB spec_color;

	base_mat_t(int tid_=-1, colorRGBA const &color_=ALPHA0, colorRGBA const &sc=BLACK, float shine_=0.0)
		: tid(tid_), shine(shine_), color(color_), spec_color(sc) {}
	bool operator==(base_mat_t const &m) const {
		return (tid == m.tid && shine == m.shine && color == m.color && spec_color == m.spec_color);
	}
};


unsigned char const SWAP_TCS_XY    = 0x01; // swap texture x and y
unsigned char const SWAP_TCS_NM_BS = 0x02; // swap normal map bitangent sign

struct obj_layer : public base_mat_t { // size = 76

	bool draw, is_emissive;
	unsigned char swap_tcs, cobj_type;
	float elastic, tscale, tdx, tdy, refract_ix, light_atten, density; // Note: elastic is misnamed - it's really hardness
	int normal_map;
	collision_func coll_func;

	obj_layer(float e=0.0, colorRGBA const &c=WHITE, bool d=0, const collision_func cf=NULL, int ti=-1, float ts=1.0, float spec=0.0, float shi=0.0)
		: base_mat_t(ti, c, colorRGB(spec, spec, spec), shi), draw(d), is_emissive(0), swap_tcs(0), cobj_type(COBJ_TYPE_STD),
		elastic(e), tscale(ts), tdx(0.0), tdy(0.0), refract_ix(1.0), light_atten(0.0), density(1.0), normal_map(-1), coll_func(cf) {}

	// assumes obj_layer contained classes are POD with no padding
	bool operator==(obj_layer const &layer) const {return (memcmp(this, &layer, sizeof(obj_layer)) == 0);}
	bool operator< (obj_layer const &layer) const {return (memcmp(this, &layer, sizeof(obj_layer)) <  0);}
	bool operator!=(obj_layer const &layer) const {return !operator==(layer);}
	bool no_draw()        const {return (!draw || color.alpha == 0.0);}
	bool has_alpha_texture() const;
	bool is_semi_trans()  const {return (color.alpha < 1.0 || has_alpha_texture());}
	bool might_be_drawn() const {return (draw || cobj_type != COBJ_TYPE_STD);}
	bool is_glass(bool shatterable=0) const {return ((tid < 0 || tid == WHITE_TEX) && (color.alpha <= 0.5 || shatterable));}
	bool swap_txy()       const {return ((swap_tcs & SWAP_TCS_XY   ) != 0);}
	bool negate_nm_bns()  const {return ((swap_tcs & SWAP_TCS_NM_BS) != 0);}
	void set_swap_tcs_flag(unsigned char mask, bool val=1) {if (val) {swap_tcs |= mask;} else {swap_tcs &= (~mask);}}
};


unsigned const COBJ_DYNAMIC     = 0x01;
unsigned const COBJ_DESTROYABLE = 0x02;
unsigned const COBJ_NO_COLL     = 0x04;
unsigned const COBJ_MOVABLE     = 0x08;
unsigned const COBJ_WAS_CUBE    = 0x10;
unsigned const COBJ_IS_INDOORS  = 0x20;

struct cobj_params : public obj_layer { // size = 84

	int cf_index;
	unsigned char surfs, flags, destroy_prob;
	//obj_layer *layer;

	cobj_params() : cf_index(-1), surfs(0), flags(0), destroy_prob(0) {}
	cobj_params(float e, colorRGBA const &c, bool d, bool id, const collision_func cf=NULL, int ci=0, int ti=-1, float ts=1.0, int s=0,
		float spec=0.0, float shi=0.0, bool nc=0) : obj_layer(e, c, d, cf, ti, ts, spec, shi), cf_index(ci), surfs(s), flags(0), destroy_prob(0) {
			if (id) {flags |= COBJ_DYNAMIC;}
			if (nc) {flags |= COBJ_NO_COLL;}
		}
};


class cobj_draw_buffer {

	obj_layer last_layer;
	vector<vert_norm_texp> verts;
	vector<vert_norm_tc> tc_verts;

public:
	int is_wet; // 0=no, 1=yes, 2=unknown

	cobj_draw_buffer() : is_wet(2) {} // initial value of is_wet is unknown
	void add_vert(vert_norm_texp const &vnt) {verts.push_back(vnt);}
	void add_vert(vert_norm const &vn, texgen_params_t const &tp) {verts.push_back(vert_norm_texp(vn, tp));}
	void add_vert(vert_norm_tc const &vntc) {tc_verts.push_back(vntc);}
	bool on_new_obj_layer(obj_layer const &l);
	void clear() {verts.clear(); tc_verts.clear(); is_wet = 2;}
	void draw() const {draw_verts(verts, GL_TRIANGLES); draw_verts(tc_verts, GL_TRIANGLES);}
	void flush() {draw(); clear();}
};


class coll_obj_group;
class csg_cube;


class coll_obj : public cube_t { // size = 244

public:
	char type, destroy, status;
	unsigned char last_coll, coll_type;
	bool fixed, is_billboard, falling;
	cobj_params cp; // could store unique cps in a set of material properties to reduce memory requirements slightly
	float radius, radius2, thickness, volume, v_fall;
	int counter, id;
	short platform_id, group_id, cgroup_id, dgroup_id, waypt_id, npoints;
	point points[N_COLL_POLY_PTS];
	vector3d norm, texture_offset;
	vector<int> occluders;

	coll_obj() : type(COLL_NULL), destroy(NON_DEST), status(COLL_UNUSED), last_coll(0), coll_type(0), fixed(0), is_billboard(0),
		falling(0), radius(0.0), radius2(0.0), thickness(0.0), volume(0.0), v_fall(0.0), counter(0), id(-1), platform_id(-1),
		group_id(-1), cgroup_id(-1), dgroup_id(-1), waypt_id(-1), npoints(0), norm(zero_vector), texture_offset(zero_vector) {}
	void init();
	void clear_internal_data();
	void setup_internal_state();
	void calc_volume();
	void calc_bcube();
	float calc_min_dim() const;
	bool clip_in_2d(float const bb[2][2], float &ztop, int d1, int d2, int dir) const;
	void set_npoints();
	void set_from_pts(point const *const pts, unsigned npts);
	void print_bounds() const;
	void bb_union(float bb[3][2], int init);
	bool is_cobj_visible() const;
	bool check_pdu_visible(pos_dir_up const &pdu) const;
	void setup_cobj_sc_texgen(vector3d const &dir, shader_t &shader) const;
	void draw_cobj(unsigned &cix, int &last_tid, int &last_group_id, shader_t &shader, cobj_draw_buffer &cdb) const;
	void get_shadow_triangle_verts(vector<vert_wrap_t> &verts, int ndiv, bool skip_spheres=0) const;
	void add_to_vector(coll_obj_group &cobjs, int type_);
	void check_if_cube();
	void add_as_fixed_cobj();
	int  add_coll_cobj();
	void re_add_coll_cobj(int index, int remove_old=1);
	bool subtract_from_cobj(coll_obj_group &new_cobjs, csg_cube const &cube, bool include_polys);
	int  intersects_cobj(coll_obj const &c, float toler=0.0) const;
	void get_side_polygons(vector<tquad_t> &sides, int top_bot_only=0) const;
	void get_contact_points(coll_obj const &c, vector<point> &contact_pts, bool vert_only=0, float toler=0.0) const;
	int  is_anchored() const;
	void move_cobj(vector3d const &vd, bool update_colls=1);
	void shift_by(vector3d const &vd, bool force=0, bool no_texture_offset=0);
	void rotate_about(point const &pt, vector3d const &axis, float angle, bool do_re_add=1);
	void add_to_platform() const;
	bool cobj_plane_side_test(point const *pts, unsigned npts, point const &lpos) const;
	bool operator<(const coll_obj &cobj) const {return (volume < cobj.volume);} // sort by size
	bool equal_params(const coll_obj &c) const {return (type == c.type && status == c.status && platform_id == c.platform_id && group_id == c.group_id && cp == c.cp);}
	bool no_draw()        const {return (status == COLL_UNUSED || status == COLL_FREED || cp.no_draw());}
	bool disabled()       const {return (status != COLL_DYNAMIC && status != COLL_STATIC);}
	bool no_collision()   const {return (disabled() || (cp.flags & COBJ_NO_COLL));}
	bool is_semi_trans()  const {return cp.is_semi_trans();}
	bool freed_unused()   const {return (status == COLL_FREED || status == COLL_UNUSED);}
	bool is_occluder()    const;// {return (status == COLL_STATIC && type == COLL_CUBE && cp.draw && !is_semi_trans());}
	bool is_big_occluder()const {return (is_occluder() && fixed && (type == COLL_POLYGON || volume > 0.001));}
	bool maybe_is_moving()const {return (platform_id >= 0 || falling);}
	bool is_movable()     const {return ((cp.flags & COBJ_MOVABLE) != 0);}
	bool is_wet()         const {return (!(cp.flags & COBJ_IS_INDOORS) && is_ground_wet());}
	bool is_snow_cov()    const {return (!(cp.flags & COBJ_IS_INDOORS) && is_ground_snowy());}
	bool may_be_dynamic() const {return (status != COLL_STATIC || maybe_is_moving() || is_movable());}
	bool is_player()      const;
	bool is_invis_player()const;
	bool truly_static()   const;
	bool no_shadow_map()  const {return (no_draw() || status != COLL_STATIC || cp.color.alpha < MIN_SHADOW_ALPHA || maybe_is_moving());}
	bool is_cylinder()    const {return (type == COLL_CYLINDER || type == COLL_CYLINDER_ROT);}
	bool is_thin_poly()   const {return (type == COLL_POLYGON && thickness <= MIN_POLY_THICK);}
	bool is_tree_leaf()   const {return is_billboard;} // we assume that a billboard cobj is a tree leaf
	bool is_cylin_vertical() const {return (points[0].x == points[1].x && points[0].y == points[1].y);}
	bool has_z_normal()      const {return (norm.x == 0.0 && norm.y == 0.0);}
	bool has_hard_edges()    const {return (type == COLL_CUBE || type == COLL_POLYGON);}
	bool has_flat_top_bot()  const {return (type == COLL_CUBE || type == COLL_POLYGON || type == COLL_CYLINDER);}
	// allow destroyable and transparent objects, drawn or opaque model3d shapes
	bool can_be_scorched()const {return (status == COLL_STATIC && !cp.has_alpha_texture() && (!no_draw() || (cp.cobj_type != COBJ_TYPE_STD && cp.color.A == 1.0)) && dgroup_id < 0);}
	point get_center_pt() const;
	point get_center_of_mass(bool ignore_group=0) const;
	float get_max_dim()   const;
	float get_light_transmit(point v1, point v2) const;
	float get_mass()      const {return volume*cp.density;}
	float get_group_mass()const;
	colorRGBA get_avg_color() const {return ((cp.tid >= 0) ? cp.color.modulate_with(texture_color(cp.tid)) : cp.color);}
	void bounding_sphere(point &center, float &brad) const;
	cylinder_3dw get_bounding_cylinder() const;
	bool is_billboard_cobj() const {return (is_billboard && is_thin_poly() && npoints == 4);}
	bool has_poly_billboard_alpha() const {return (is_billboard_cobj() && cp.has_alpha_texture());}
	bool check_poly_billboard_alpha(point const &p1, point const &p2, float t) const;
	bool line_intersect(point const &p1, point const &p2) const;
	bool line_int_exact(point const &p1, point const &p2, float &t, vector3d &cnorm, float tmin=0.0, float tmax=1.0) const;
	bool sphere_intersects_exact(point const &sc, float sr, vector3d &cnorm, point &new_sc) const;
	bool intersects_all_pts(point const &pos, point const *const pts, unsigned npts) const; // coll_cell_search.cpp
	void convert_cube_to_ext_polygon();
	colorRGBA get_color_at_point(point const &pos, vector3d const &normal, bool fast) const;
	bool is_occluded_from_camera() const;
	void register_coll(unsigned char coll_time, unsigned char coll_type_);
	void create_portal() const; // destroy_cobj.cpp
	void add_connect_waypoint(); // waypoints.cpp
	void remove_waypoint();
	void check_indoors_outdoors();
	void write_to_cobj_file(std::ostream &out, coll_obj &prev) const;
	
	// inexact primitive intersections
	int cube_intersects(cube_t const &cube) const;
	int sphere_intersects(point const &pos, float radius) const;
	int sphere_intersects(sphere_t const &sphere) const {return sphere_intersects(sphere.pos, sphere.radius);}
	int contains_point(point const &pos) const;

	vector3d get_cobj_supporting_normal(point const &support_pos, bool bot_surf=0) const;
	vector3d get_cobj_resting_normal() const;
	bool is_point_supported(point const &pos) const;

	// drawing code
	void setup_cube_face_texgen(texgen_params_t &tp, unsigned tdim0, unsigned tdim1, float const tscale[2]) const;
	void draw_coll_cube(int tid, cobj_draw_buffer &cdb) const;
	void set_poly_texgen(int tid, vector3d const &normal, shader_t &shader) const;
	void get_polygon_tparams(int tid, vector3d const &normal, texgen_params_t &tp) const;
	void draw_polygon(int tid, point const *pts, int npts, vector3d const &normal, cobj_draw_buffer &cdb) const;
	void draw_extruded_polygon(int tid, cobj_draw_buffer &cdb) const;
	void draw_cylin_ends(int tid, int ndiv, cobj_draw_buffer &cdb) const;
};


struct cobj_id_set_t : public set<unsigned> {

	void must_insert(unsigned index) {
		bool const did_ins(insert(index).second);
		assert(did_ins);
	}
	void must_erase(unsigned index) {
		unsigned const num_removed(erase(index));
		assert(num_removed == 1);
	}
	bool has(unsigned index) {return find(index) != end();}
};

struct cgroup_props_t {

	bool valid;
	float volume, mass;
	point center_of_mass;
	cgroup_props_t() : valid(0), volume(0.0), mass(0.0), center_of_mass(all_zeros) {}
};

struct cobj_group_t : public cobj_id_set_t, public cgroup_props_t {
	
	// Note inserting and erasing invalidate the props, forcing recalculation any time someone calls get_props()
	void must_insert(unsigned index) {cobj_id_set_t::must_insert(index); valid = 0;}
	void must_erase (unsigned index) {cobj_id_set_t::must_erase (index); valid = 0;}
	void update_props();
	cgroup_props_t const &get_props() {if (!valid) {update_props();} return *this;}
};

struct cobj_groups_t : public vector<cobj_group_t> {

	unsigned new_group() {unsigned const id(size()); push_back(cobj_group_t()); return id;}
	void invalidate_group(unsigned gid)               {assert(gid < size()); operator[](gid).valid = 0;}
	void add_cobj   (unsigned gid, unsigned cid)      {assert(gid < size()); operator[](gid).must_insert(cid);}
	void remove_cobj(unsigned gid, unsigned cid)      {assert(gid < size()); operator[](gid).must_erase (cid);}
	cobj_id_set_t  const &get_set(unsigned gid) const {assert(gid < size()); return operator[](gid);}
	cgroup_props_t const &get_props(unsigned gid)     {assert(gid < size()); return operator[](gid).get_props();}
};


class coll_obj_group : public vector<coll_obj> {

public:
	bool has_lt_atten, has_voxel_cobjs;
	cobj_id_set_t dynamic_ids, drawn_ids, platform_ids;

	coll_obj_group() : has_lt_atten(0), has_voxel_cobjs(0) {}
	void clear_ids();
	void clear();
	void finalize();
	void remove_invalid_cobjs();
	void check_cubes();
	void merge_cubes();
	void process_negative_shapes();
	void remove_overlapping_cubes(int min_split_destroy_thresh);
	void subdiv_cubes();
	void sort_cobjs_for_rendering();
	void set_coll_obj_props(int index, int type, float radius, float radius2, int platform_id, cobj_params const &cparams);
	void remove_index_from_ids(int index);
	
	coll_obj &get_cobj(int index) {
		assert(index >= 0 && index < (int)size());
		return operator[](index);
	}
};


class cobj_draw_groups {

	struct cobj_draw_group {
		int parent; // -1 is unset
		point parent_pos;
		vector<unsigned> ids; // typically a contiguous range within dcobjs

		cobj_draw_group() : parent(-1), parent_pos(all_zeros) {}
		bool empty() const {return ids.empty();}
	};

	vector<coll_obj> dcobjs;
	vector<cobj_draw_group> groups;

	cobj_draw_group       &get_group(int group_id)       {assert(group_id >= 0 && group_id < (int)groups.size()); return groups[group_id];}
	cobj_draw_group const &get_group(int group_id) const {assert(group_id >= 0 && group_id < (int)groups.size()); return groups[group_id];}

public:
	coll_obj const &get_cobj(unsigned index) const {
		assert(index < dcobjs.size());
		return dcobjs[index];
	}
	unsigned new_group();
	void add_to_group(coll_obj const &cobj); // Note: cobj.dgroup_id should be set correctly
	bool set_parent_or_add_cobj(coll_obj const &cobj); // Note: cobj.dgroup_id should be set correctly
	vector<unsigned> const &get_draw_group(int group_id, coll_obj const &parent);
};


struct cobj_query_callback {
	virtual bool register_cobj(coll_obj const &cobj) = 0;
};


class polygon_t : public vector<vert_norm_tc> {

public:
	colorRGBA color;

	polygon_t(colorRGBA const &c=ALPHA0) : color(c) {}
	polygon_t(vector<vert_norm_tc> const &vv, colorRGBA const &c=ALPHA0) : vector<vert_norm_tc>(vv), color(c) {}
	polygon_t(triangle const &t, colorRGBA const &c=WHITE) : color(c) {from_triangle(t);}
	void from_triangle(triangle const &t);
	bool is_convex() const;
	bool is_coplanar(float thresh) const;
	vector3d get_planar_normal() const;
	bool is_valid() const {return (size() >= 3 && is_triangle_valid((*this)[0].v, (*this)[1].v, (*this)[2].v));}
	void from_points(vector<point> const &pts);
};


struct coll_tquad : public tquad_t { // size = 68

	vector3d normal;
	
	union {
		unsigned cid;
		color_wrapper color;
	};
	coll_tquad() {}
	coll_tquad(coll_obj const &c);
	coll_tquad(polygon_t const &p);
	coll_tquad(triangle const &t, colorRGBA const &c=WHITE);
	void update_normal() {get_normal(pts[0], pts[1], pts[2], normal, 1);}

	static bool is_cobj_valid(coll_obj const &c) {
		return (!c.disabled() && c.type == COLL_POLYGON && (c.npoints == 3 || c.npoints == 4) && c.thickness <= MIN_POLY_THICK);
	}
	bool line_intersect(point const &p1, point const &p2) const {
		float t;
		return line_poly_intersect(p1, p2, pts, npts, normal, t);
	}
	bool line_int_exact(point const &p1, point const &p2, float &t, vector3d &cnorm, float tmin, float tmax) const {
		if (!line_poly_intersect(p1, p2, pts, npts, normal, t) || t > tmax || t < tmin) return 0;
		cnorm = get_poly_dir_norm(normal, p1, (p2 - p1), t);
		return 1;
	}
};


void copy_polygon_to_cobj(polygon_t const &poly, coll_obj &cobj);
void copy_tquad_to_cobj(coll_tquad const &tquad, coll_obj &cobj);


struct coll_cell { // size = 52

	float zmin, zmax;
	vector<int> cvals;

	void clear(bool clear_vectors);

	void update_zmm(float zmin_, float zmax_) {
		assert(zmin_ <= zmax_);
		zmin = min(zmin_, zmin);
		zmax = max(zmax_, zmax);
	}
	void add_entry(int index) {
		if (INIT_CCELL_SIZE > 0 && cvals.capacity() == 0) {cvals.reserve(INIT_CCELL_SIZE);}
		cvals.push_back(index);
	}
};


struct color_tid_vol : public cube_t {

	int cid, tid, destroy;
	bool draw, unanchored, is_2d;
	float volume, thickness, tscale, max_frag_sz;
	colorRGBA color;

	color_tid_vol(coll_obj const &cobj, float volume_, float thickness_, bool ua);
	bool maybe_is_glass() const;
};


class platform { // animated (player controlled) scene object

	// constants
	bool const cont; // continuous - always in motion
	bool const is_rot; // motion is rotation rather than translation
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

	multi_trigger_t triggers;

	void move_platform(float dist_traveled);
	bool empty() const {return (cobjs.empty() && lights.empty());}

public:
	// other data
	vector<unsigned> cobjs;  // collision object(s) bound to this platform
	vector<unsigned> lights; // dynamic light source(s) bound to this platform
	
	platform(float fs=1.0, float rs=1.0, float sd=0.0, float rd=0.0, float dst=1.0, float ad=0.0,
		point const &o=all_zeros, vector3d const &dir_=plus_z, bool c=0, bool ir=0);
	void add_triggers(multi_trigger_t const &t) {triggers.add_triggers(t);} // deep copy
	bool has_dynamic_shadows() const {return (cont || state >= ST_FWD);}
	vector3d get_delta()       const {return (pos - origin);}
	vector3d get_range()       const {return dir*ext_dist;}
	vector3d get_last_delta()  const {return delta;}
	vector3d get_velocity()    const;
	void clear_last_delta() {delta = all_zeros;}
	void add_cobj(unsigned cobj);
	void add_light(unsigned light) {lights.push_back(light);}
	void next_frame();
	void shift_by(vector3d const &val);
	void reset();
	void activate();
	bool check_activate(point const &p, float radius, int activator);
	void advance_timestep();
	bool is_moving  () const {return (state == ST_FWD || state == ST_REV);}
	bool is_active  () const {return (state != ST_NOACT);}
	bool is_rotation() const {return is_rot;}
	void write_to_cobj_file(std::ostream &out) const;
};


struct platform_cont : public deque<platform> {

	bool add_from_file(FILE *fp, geom_xform_t const &xf, multi_trigger_t const &triggers);
	void check_activate(point const &p, float radius, int activator);
	void shift_by(vector3d const &val);
	void add_current_cobjs();
	void advance_timestep();
	bool any_active() const;
	bool any_moving_platforms_in_view(pos_dir_up const &pdu) const;
	
	platform &get_cobj_platform(coll_obj const &cobj) {
		assert(cobj.platform_id >= 0 && cobj.platform_id < (int)size());
		return operator[](cobj.platform_id);
	}
};


struct shadow_sphere : public sphere_t {

	bool is_player;
	char ctype;
	int cid;

	shadow_sphere() : is_player(0), ctype(COLL_NULL), cid(-1) {}
	shadow_sphere(point const &pos0, float radius0, int cid0, bool is_player_=0);
	
	inline bool line_intersect(point const &p1, point const &p2) const {
		return (line_sphere_intersect(p1, p2, pos, radius) && (ctype == COLL_SPHERE || line_intersect_cobj(p1, p2)));
	}
	bool line_intersect_cobj(point const &p1, point const &p2) const;

	inline bool test_volume(point const *const pts, unsigned npts, point const &lpos) const {
		return (cid < 0 || test_volume_cobj(pts, npts, lpos));
	}
	bool test_volume_cobj(point const *const pts, unsigned npts, point const &lpos) const;
};


class obj_draw_group { // vbo_wrap_t?

	unsigned start_cix, end_cix, vbo, num_verts;
	bool use_vbo, inside_beg_end;
	vector<vert_norm> verts;

public:
	obj_draw_group(bool use_vbo_=0) : start_cix(0), end_cix(0), vbo(0), num_verts(0), use_vbo(use_vbo_), inside_beg_end(0) {}
	void free_vbo();
	void draw_vbo() const;
	bool begin_render(unsigned &cix);
	void end_render();
	void add_draw_polygon(point const *const points, vector3d const &normal, unsigned npoints, unsigned cix);
	void set_vbo_enable(bool en) {assert(vbo == 0); use_vbo = en;} // can't call while the vbo is valid
	bool vbo_enabled() const {return use_vbo;}
};


#endif // _COLLISION_DETECT_H_

