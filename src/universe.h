// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/4/02

#ifndef _UNIVERSE_H_
#define _UNIVERSE_H_

#include "3DWorld.h"
#include "universe_base.h"
#include "upsurface.h"
#include "draw_utils.h"
#include <map>
#include <sstream>

using std::string;
using std::ostringstream;
using std::ostream;
using std::istream;


class s_object;
class uasteroid;
class uasteroid_field;
class uasteroid_belt;
class uasteroid_belt_system;
class uasteroid_belt_planet;


// stellar object types - must be ordered largest to smallest
enum {UTYPE_NONE=0, UTYPE_CELL, UTYPE_GALAXY, UTYPE_SYSTEM, UTYPE_STAR, UTYPE_PLANET, UTYPE_MOON, UTYPE_SURFACE, UTYPE_ASTEROID, NUM_UTYPES};

float const u_exp_size[NUM_UTYPES] = {
	  0.0,          0.0,        0.0,          0.0,          6.0,        4.0,          3.0,        0.0,           2.0
};

// stellar object modifications
enum {MOD_DESTROYED=0, MOD_OWNER, MOD_NAME, N_UMODS};

int const AST_BELT_ID = -2;


float const GALAXY_SCALE    = 8.0;
float const GALAXY_OVERLAP  = 0.5;
float const GALAXY_MIN_SIZE = 18.0*GALAXY_SCALE;
float const GALAXY_MAX_SIZE = 24.0*GALAXY_SCALE;
float const STAR_MAX_SIZE   = 0.14;
float const STAR_MIN_SIZE   = 0.04;
float const PLANET_MAX_SIZE = 0.035;
float const PLANET_MIN_SIZE = 0.008;
float const MOON_MAX_SIZE   = 0.008;
float const MOON_MIN_SIZE   = 0.003;

float const SYSTEM_MIN_SPACING         = 5.0;
float const PLANET_TO_SUN_MIN_SPACING  = 0.25;
float const PLANET_TO_SUN_MAX_SPACING  = 1.8;
float const INTER_PLANET_MIN_SPACING   = 0.04;
float const MOON_TO_PLANET_MIN_SPACING = 0.025;
float const MOON_TO_PLANET_MAX_SPACING = 0.15;
float const MOON_TO_PLANET_MIN_GAP     = 0.008;
float const INTER_MOON_MIN_SPACING     = 0.01;
float const MIN_RAD_SPACE_FACTOR       = 1.2;

float const MP_COLOR_VAR        = 0.4;
float const PLANET_ATM_RSCALE   = 1.025;
float const ORBIT_PLANE_DELTA   = 0.06;
float const ORBIT_SPACE_MARGIN  = 1.1;
float const RANGE_OFFSET        = 0.0001;
float const MOON_HMAP_SCALE     = 0.08;
float const PLANET_HMAP_SCALE   = 0.04;
float const SURFACE_HEIGHT      = 0.01; // multiples of the radius
float const INITIAL_FREQ        = 0.7;
float const TEX_HEIGHT_SCALE    = 1.0;
float const NEAR_CLIP_SCALE     = 1.5;

float const BASE_AMBIENT  = 0.25;
float const BASE_DIFFUSE  = 1.3;
float const WHITE_COMP_A  = 0.5;
float const ATTEN_AMB_VAL = 0.7;
float const WHITE_COMP_D  = 0.1;
float const STAR_CONST    = 0.4;
float const STAR_LINEAR   = 2.4;
float const STAR_QUAD     = 0.0;

float const OM_WCA(1.0 - WHITE_COMP_A);
float const OM_WCD(1.0 - WHITE_COMP_D);
float const OM_AAV(1.0 - ATTEN_AMB_VAL);
float const STAR_LINEAR_SCALED(STAR_LINEAR/PLANET_TO_SUN_MAX_SPACING);
float const STAR_QUAD_SCALED(STAR_QUAD/(PLANET_TO_SUN_MAX_SPACING*PLANET_TO_SUN_MAX_SPACING));
float const MAX_PLANET_EXTENT(MOON_TO_PLANET_MAX_SPACING + MOON_MAX_SIZE);
float const MAX_SYSTEM_EXTENT(PLANET_TO_SUN_MAX_SPACING + MAX_PLANET_EXTENT);
float const MAX_GALAXY_EXTENT(GALAXY_MAX_SIZE + MAX_SYSTEM_EXTENT);
float const AVG_STAR_SIZE(0.5*(STAR_MAX_SIZE + STAR_MIN_SIZE));

unsigned const U_BLOCKS      = 7;
unsigned const U_BLOCKS_SQ   = U_BLOCKS*U_BLOCKS;
unsigned const U_BLOCKS_CU   = U_BLOCKS_SQ*U_BLOCKS;
unsigned const U_BLOCKSo2    = (U_BLOCKS-1)/2;

float const CELL_EDGE        = CELL_SIZE*((U_BLOCKS-1.0)/2.0 + 0.5);
float const CELL_SIZEo2      = CELL_SIZE/2.0;
float const CELL_SPHERE_RAD  = CELL_SIZEo2*sqrt(3.0);
float const U_VIEW_DIST      = U_BLOCKS*CELL_SIZEo2;
float const CELL_SIZE_INV    = 1.0/CELL_SIZE;
float const RS_SCALE         = 7.0*CELL_SIZE_INV;
float const TEX_H_SCALE      = TEX_HEIGHT_SCALE/(SURFACE_HEIGHT*sqrt((double)SINES_PER_FREQ));
float const NEAR_CLIP_SCALED = NEAR_CLIP_SCALE*UNIV_NCLIP_SCALE*NEAR_CLIP;

colorRGBA const P_WATER_C(0.2, 0.2, 0.8, 1.0);
colorRGBA const P_ICE_C(  0.5, 0.7, 0.9, 1.0);


// forward references
class umoon;
class uplanet;
class ustar;
class ussystem;
class ugalaxy;
class ucell;
class ushader_group;

struct cell_ixs_t {int ix[3];};


struct shadow_vars_t {
	point sun_pos, ss_pos;
	vector3d rscale;
	float sun_radius, ss_radius, ring_ri, ring_ro;

	shadow_vars_t() : sun_radius(0.0), ss_radius(0.0), ring_ri(0.0), ring_ro(0.0) {}
	shadow_vars_t(point const &sp, float sr, point const &ssp, float ssr, vector3d const &rs, float rri, float rro)
		: sun_pos(sp), ss_pos(ssp), sun_radius(sr), ss_radius(ssr), rscale(rs), ring_ri(rri), ring_ro(rro) {}
};


class named_obj { // size = 16

	string name;

public:
	named_obj() : name("Unnamed") {}
	named_obj(string const &name_) : name(name_) {}
	void setname(string const &name_) {name = name_;}
	string const &getname() const {return name;}
	void gen_name(s_object const &sobj);
	bool rename(s_object const &sobj, string const &name_);
	bool lookup_given_name(s_object const &sobj);
};


class uobj_rgen: public uobject { // size = 56

public:
	char gen;
	int urseed1, urseed2;

	uobj_rgen() : gen(0), urseed1(1), urseed2(1) {}
	void gen_rseeds();
	void get_rseeds();
	void set_rseeds() const;
	int get_id() const {return urseed1;} // not complete id, but should be good enough
};


class uobj_solid : public uobj_rgen, public named_obj { // size = 176

public:
	char type;
	unsigned short num_satellites;
	float temp, density, mass, gravity;
	colorRGBA color, colorA, colorB;

	uobj_solid(char type_) : type(type_), num_satellites(0), temp(0.0), density(0.0), mass(0.0), gravity(0.0) {set_defaults();}
	virtual ~uobj_solid() {}
	void set_defaults() {status = 0; gen = 0;}
	void get_colors(unsigned char ca[3], unsigned char cb[3]) const;
	void adjust_colorAB(float delta);
	void gen_colorAB(float delta);
	void set_grav_mass();
	bool collision(point const &p, float rad, vector3d const &v, point &cpos, float &coll_r, bool simple) const;
	bool rename(std::string const &name_) {setname(name_); return 1;}
	int  get_type() const {return type;}

	void add_gravity_vector(vector3d &vgravity, point const &mpos) const { // inlined
		add_gravity_vector_base(vgravity, mpos, gravity, MAX_SOBJ_GRAVITY);
	}
	virtual bool surface_test(float rad, point const &p, float &coll_r, bool simple) const {return 1;}
	virtual void free_uobj() {gen = 0;}
};


struct rotated_obj { // size = 20

	vector3d rot_axis;
	double rev_ang, rev_ang0, rot_ang, rot_ang0; // Note: rev_ang here to preserve rand() call order

	rotated_obj() : rev_ang(0.0), rev_ang0(0.0), rot_ang(0.0), rot_ang0(0.0) {}
	void rgen_values();
	void apply_gl_rotate() const;
	void rotate_vector(vector3d &v) const;
	void rotate_vector_inv(vector3d &v) const;
};


class urev_body : public uobj_solid, public color_gen_class, public rotated_obj { // size = 268

	// for textures/colors
	unsigned char a[3], b[3];
	
protected:
	void calc_snow_thresh();

public:
	bool gas_giant; // planets only?
	int owner;
	unsigned orbiting_refs, tid, tsize;
	float orbit, rot_rate, rev_rate, atmos, water, lava, resources, cloud_scale, wr_scale, snow_thresh;
	vector3d rev_axis, v_orbit;
	std::shared_ptr<upsurface> surface;
	string comment;

	urev_body(char type_) : uobj_solid(type_), gas_giant(0), owner(NO_OWNER), orbiting_refs(0), tid(0), tsize(0), orbit(0.0), rot_rate(0.0),
		rev_rate(0.0), atmos(0.0), water(0.0), lava(0.0), resources(0.0), cloud_scale(1.0), wr_scale(1.0), snow_thresh(0.0) {}
	virtual ~urev_body() {unset_owner();}
	void gen_rotrev();
	template<typename T> bool create_orbit(vector<T> const &objs, int i, point const &pos0, vector3d const &raxis,
		float radius0, float max_size, float min_size, float rspacing, float ispacing, float minspacing, float min_gap);
	void gen_surface();
	void check_gen_texture(unsigned size);
	void create_rocky_texture(unsigned size);
	void create_gas_giant_texture();
	void gen_texture_data_and_heightmap(unsigned char *data, unsigned size);
	bool has_heightmap() const {return (surface != nullptr && surface->has_heightmap() && !use_procedural_shader());}
	bool surface_test(float rad, point const &p, float &coll_r, bool simple) const;
	float get_radius_at(point const &p, bool exact=0) const;
	float get_dheight_at(point const &p, bool exact=0) const;
	virtual bool has_custom_shadow_profile() const {return (surface != nullptr);}
	bool pt_over_land(point const &p) const;
	bool can_land() const;
	bool can_land_at(point const &p) const {return (can_land() && pt_over_land(p));}
	bool colonizable() const;
	bool liveable() const;
	float get_land_value(unsigned align, point const &cur_pos, float sradius) const;
	string get_info() const;
	void set_owner(s_object const &sobj, unsigned owner_);
	void set_owner_int(unsigned owner_);
	void unset_owner();
	void check_owner(s_object const &sobj);
	void get_owner_info(ostringstream &oss) const;
	int  get_owner() const {return owner;}
	colorRGBA get_owner_color() const;
	bool use_procedural_shader() const;
	void upload_colors_to_shader(shader_t &s) const;
	void get_surface_color(unsigned char *data, float val, float phi) const;
	bool draw(point_d pos_, ushader_group &usg, pt_line_drawer planet_plds[2], shadow_vars_t const &svars, bool use_light2);
	void draw_surface(point_d const &pos_, float radius0, float size, int ndiv);
	void show_colonizable_liveable(point const &pos_, float radius0) const;
	void inc_orbiting_refs() {++orbiting_refs;}
	void dec_orbiting_refs(s_object const &sobj);
	void explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
		int align, unsigned eflags=0, free_obj const *parent_=NULL);
	
	virtual float get_hmap_scale() const = 0;
	virtual void create(bool phase) = 0;
	virtual void calc_temperature() = 0;
	virtual void get_valid_orbit_r(float &orbit_r, float obj_r) const = 0;
	virtual bool colonizable_int() const = 0;
	virtual float get_vegetation() const = 0;
	virtual point_d do_update(point_d const &p0, bool update_rev=1, bool update_rot=1);
	virtual void free_texture();
	virtual void free_uobj();
};


class uplanet : public urev_body { // size = 324
public:
	unsigned population; // unused
	float mosize, ring_ri, ring_ro;
	vector3d rscale;
	vector<umoon> moons;
	vector<color_wrapper> ring_data;
	ussystem *system;
	std::shared_ptr<uasteroid_belt_planet> asteroid_belt;
	unsigned ring_tid;
	// trade items?

	uplanet() : urev_body(UTYPE_PLANET), population(0), mosize(0.0), ring_ri(0.0), ring_ro(0.0), system(NULL), ring_tid(0) {}
	void create(bool phase);
	void process();
	point_d do_update(point_d const &p0, bool update_rev=1, bool update_rot=1);
	void gen_prings();
	void gen_color();
	void calc_temperature();
	void get_valid_orbit_r(float &orbit_r, float obj_r) const;
	bool colonizable_int() const {return (!gas_giant && water > 0.1 && atmos > 0.1);}
	bool has_vegetation()  const {return (!gas_giant && water > 0.2 && atmos > 0.1);}
	float get_vegetation() const;
	void ensure_rings_texture();
	void draw_prings(ushader_group &usg, upos_point_type const &pos_, float size_, point const &sun_pos, float sun_radius, bool dir) const;
	void draw_atmosphere(ushader_group &usg, upos_point_type const &pos_, float size_, shadow_vars_t const &svars) const;
	void free_texture();
	void free_uobj();
	float get_hmap_scale () const {return PLANET_HMAP_SCALE;}
	float get_ring_rscale() const {return max(rscale.x, rscale.y)*ring_ro/radius;}
	string get_name() const {return "Planet " + getname();}
};


class umoon : public urev_body { // size = 268

public:
	uplanet *planet;

	umoon() : urev_body(UTYPE_MOON), planet(NULL) {}
	void create(bool phase);
	void gen_color();
	void calc_temperature();
	bool shadowed_by_planet();
	void get_valid_orbit_r(float &orbit_r, float obj_r) const;
	bool colonizable_int() const {return (radius > 1.5*MOON_MIN_SIZE && planet && planet->colonizable());}
	float get_vegetation() const {return 0;}
	float get_hmap_scale() const {return MOON_HMAP_SCALE;}
	string get_name() const {return "Moon " + getname();}
};


class ustar : public uobj_solid { // size = 176

	struct solar_flare {
		float length, radius, angle;
		unsigned time, lifetime;
		vector3d dir;
		colorRGBA color1, color2;

		solar_flare() : length(0.0), radius(0.0), angle(0.0), time(0), lifetime(0) {}
		void gen(colorRGBA const &color);
		void update(colorRGBA const &color);
		void draw(float size, int ndiv, bool texture) const;
	};

	vector<solar_flare> solar_flares;

public:
	vector3d rot_axis;

	ustar() : uobj_solid(UTYPE_STAR) {}
	void create(point const &pos_);
	void gen_color();
	colorRGBA get_ambient_color_val() const;
	colorRGBA get_light_color() const;
	bool draw(point_d pos_, ushader_group &usg, pt_line_drawer_no_lighting_t &star_pld, point_sprite_drawer &star_psd, bool distant);
	void draw_flares(int ndiv, bool texture);
	float get_energy() const {return (is_ok() ? PLANET_TO_SUN_MAX_SPACING*PLANET_TO_SUN_MAX_SPACING*temp*radius : 0.0);}
	vector3d get_solar_wind_accel(point const &obj_pos, float obj_mass, float obj_surf_area) const;
	string get_name()  const {return "Star " + getname();}
	string get_info()  const;
	int get_owner() const {return NO_OWNER;}
	void explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
		int align, unsigned eflags=0, free_obj const *parent_=NULL);
};


class ussystem : public uobj_rgen { // size = 268

public:
	unsigned cluster_id;
	ustar sun;
	vector<uplanet> planets;
	std::shared_ptr<uasteroid_belt_system> asteroid_belt;
	ugalaxy *galaxy;
	colorRGBA galaxy_color;
	
	ussystem() : cluster_id(0), galaxy(NULL), galaxy_color(ALPHA0) {}
	void create(point const &pos_);
	void calc_color();
	void process();
	colorRGBA const &get_galaxy_color();
	void free_uobj();
	void free_planets();
	string get_name() const {return "System " + sun.getname();}
};


class unebula : public uobject_base, public volume_part_cloud {

	colorRGBA color[3];

public:
	void gen(float range, ellipsoid_t const &bounds);
	void draw(point_d pos_, point const &camera, float max_dist, vpc_shader_t &s) const;
	void free_uobj() {points.clear();}
	bool is_valid() const {return !points.empty();}
};


class ugalaxy : public uobj_rgen, public named_obj, public ellipsoid_t { // size = 148 (164)

	mutable float lrq_rad;
	mutable point lrq_pos;

	void apply_scale_transform(point &pos_) const;
	point gen_valid_system_pos() const;

public:
	struct system_cluster {

		float radius, bounds;
		point center;
		colorRGBA color;
		vector<point> systems;
		unsigned s1, s2;

		system_cluster(float radius_, point const &center_) : radius(radius_), bounds(0.0), center(center_), color(BLACK), s1(0), s2(0) {}
	};

	vector<ussystem> sols;
	deque<system_cluster> clusters;
	vector<uasteroid_field> asteroid_fields;
	unebula nebula;
	colorRGBA color;

	ugalaxy();
	~ugalaxy();
	void calc_color();
	void calc_bounding_sphere();
	bool create(ucell const &cell, int index);
	float get_radius_at(point const &pos_, bool exact=0) const;
	bool is_close_to(ugalaxy const &g, float overlap_amount) const;
	void process(ucell const &cell);
	bool gen_system_loc(vector<point> const &placed);
	void clear_systems();
	void free_uobj();
	string get_name() const {return "Galaxy " + getname();}
};


class ucell : public uobj_rgen { // size = 84

	pt_line_drawer planet_plds[2]; // {1-pixel, 2-pixel}
	pt_line_drawer_no_lighting_t star_pld; // 1-pixel
	point_sprite_drawer star_psd; // 2-pixel
	colorRGBA last_bkg_color;
	point last_player_pos;
	unsigned last_star_cache_ix;
	bool cached_stars_valid;

	void draw_all_stars(ushader_group &usg, bool clear_star_pld);

public:
	point rel_center;
	std::shared_ptr<vector<ugalaxy> > galaxies; // must be a pointer to a vector to avoid deep copies

	ucell() : last_bkg_color(BLACK), last_player_pos(all_zeros), last_star_cache_ix(0), cached_stars_valid(0) {}
	void gen_cell(int const ii[3]);
	void draw_nebulas(ushader_group &usg) const;
	void draw_systems(ushader_group &usg, s_object const &clobj, unsigned pass, bool no_move, bool skip_closest, bool sel_cell, bool gen_only);
	void free_uobj();
	bool is_visible() const;
	string get_name() const {return "Universe Cell";}
	void free_context() {star_pld.free_vbo();}
};


class s_object { // size = 56

public:
	int cellxyz[3], galaxy, cluster, system, planet, moon, asteroid_field, asteroid, type, val, id;
	float dist;
	uobj_solid *object;

	s_object() {init();}
	bool write(ostream &out) const;
	bool read(istream &in);
	void init();
	void assign(int gc, int cl, int sy, float di, int ty, uobj_solid *obj);
	void assign_asteroid(float d, int ai, int afi) {type = UTYPE_ASTEROID; asteroid_field = afi; asteroid = ai; dist = d;}
	bool operator<(const s_object &I) const;
	bool bad_cell() const;
	bool is_solid() const {return (type == UTYPE_STAR || type == UTYPE_PLANET || type == UTYPE_MOON || type == UTYPE_ASTEROID);}
	void print() const;
	bool is_destroyed() const;
	void register_destroyed_sobj() const;
	unsigned get_owner() const;
	void set_owner(unsigned owner) const;
	bool has_valid_system() const {return (type >= UTYPE_SYSTEM && system >= 0);}
	ucell     &get_ucell()  const;

	ugalaxy   &get_galaxy() const {
		ucell &c(get_ucell());
		assert(type >= UTYPE_GALAXY && (unsigned)galaxy < c.galaxies->size());
		return (*c.galaxies)[galaxy];
	}
	ugalaxy::system_cluster &get_cluster() const {
		ugalaxy &g(get_galaxy());
		assert(type >= UTYPE_SYSTEM && (unsigned)cluster < g.clusters.size());
		return g.clusters[cluster];
	}
	ussystem &get_system() const {
		ugalaxy &g(get_galaxy());
		assert(has_valid_system() && (unsigned)system < g.sols.size());
		return g.sols[system];
	}
	ustar &get_star() const {return get_system().sun;}

	uplanet &get_planet() const {
		ussystem &s(get_system());
		assert(type >= UTYPE_PLANET && (unsigned)planet < s.planets.size());
		return s.planets[planet];
	}
	umoon &get_moon() const {
		uplanet &p(get_planet());
		assert(type >= UTYPE_MOON && (unsigned)moon < p.moons.size());
		return p.moons[moon];
	}
	urev_body &get_world() const {
		if (type == UTYPE_PLANET) return get_planet();
		assert(type == UTYPE_MOON);
		return get_moon();
	}
	uasteroid_field &get_asteroid_field() const;
	uasteroid_belt  &get_asteroid_belt () const;
	uasteroid &get_asteroid() const;
	void destroy_asteroid() const;
};


struct cell_block {

	ucell cells[U_BLOCKS][U_BLOCKS][U_BLOCKS];
};


class universe_t : protected cell_block {

public:
	void init();
	void shift_cells(int dx, int dy, int dz);
	void free_textures();
	void draw_all_cells(s_object const &clobj, bool skip_closest, bool no_move, int no_distant, bool gen_only=0);
	int get_closest_object(s_object &result, point pos, int max_level, bool include_asteroids, bool offset, float expand, bool get_destroyed=0, float g_expand=1.0) const;
	bool get_trajectory_collisions(s_object &result, point &coll, vector3d dir, point start, float dist, float line_radius) const;
	float get_point_temperature(s_object const &clobj, point const &pos, point &sun_pos) const;

	int get_object_closest_to_pos(s_object &result, point const &pos, bool include_asteroids, float expand=1.0) const {
		return get_closest_object(result, pos, UTYPE_MOON, include_asteroids, 1, expand);
	}
	int get_close_system(point const &pos, s_object &result, float expand) const {
		if (!get_closest_object(result, pos, UTYPE_SYSTEM, 0, 1, expand)) return 0; // find closest system (check last param=offset?)
		return result.has_valid_system();
	}
	bool bad_cell_xyz(int const cxyz[3]) const {
		if (cxyz[0] < 0              || cxyz[1] < 0              || cxyz[2] < 0)              return 1;
		if (cxyz[0] >= int(U_BLOCKS) || cxyz[1] >= int(U_BLOCKS) || cxyz[2] >= int(U_BLOCKS)) return 1;
		return 0;
	}
	ucell const &get_cell(int const cxyz[3]) const {
		assert(!bad_cell_xyz(cxyz));
		return cells[cxyz[2]][cxyz[1]][cxyz[0]];
	}
	ucell &get_cell(int const cxyz[3]) {
		assert(!bad_cell_xyz(cxyz));
		return cells[cxyz[2]][cxyz[1]][cxyz[0]];
	}
	ucell const &get_cell(s_object const &so) const {return get_cell(so.cellxyz);}
	ucell       &get_cell(s_object const &so)       {return get_cell(so.cellxyz);}
};


typedef string modmap_val_t;
typedef map<s_object, modmap_val_t> modmap;


struct coll_test { // size = 16

	int index;
	float dist, rad, t;
	bool operator<(const coll_test &A) const {return dist < A.dist;}
};


inline uplanet const &get_planet(s_object const &so) {return so.get_planet();}


bool import_default_modmap();
bool import_modmap(string const &filename);
bool export_modmap(string const &filename);
s_object get_shifted_sobj(s_object const &sobj);
float calc_sphere_size(point const &pos, point const &camera, float radius, float d_adj=0.0);
bool sphere_size_less_than(point const &pos, point const &camera, float radius, float num_pixels);


#endif


