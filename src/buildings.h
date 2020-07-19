// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h" // for vbo_wrap_t

bool const ADD_BUILDING_INTERIORS  = 1;
bool const EXACT_MULT_FLOOR_HEIGHT = 1;
unsigned const MAX_CYLIN_SIDES     = 36;
unsigned const MAX_DRAW_BLOCKS     = 8; // for building interiors only; currently have floor, ceiling, walls, and doors
unsigned const NUM_STAIRS_PER_FLOOR= 12;
float const FLOOR_THICK_VAL        = 0.1; // 10% of floor spacing

class light_source;
class lmap_manager_t;
class building_nav_graph_t;
struct pedestrian_t;

struct building_occlusion_state_t {
	point pos;
	vector3d xlate;
	vector<unsigned> building_ids;
	vector<point> temp_points;

	void init(point const &pos_, vector3d const &xlate_) {
		pos   = pos_;
		xlate = xlate_;
		building_ids.clear();
	}
};

struct cube_with_zval_t : public cube_t {
	float zval;
	cube_with_zval_t() : zval(0.0) {}
	cube_with_zval_t(cube_t const &c, float zval_=0.0) : cube_t(c), zval(zval_) {}
};

typedef vector<cube_with_zval_t> vect_cube_with_zval_t;

struct tid_nm_pair_t { // size=28

	int tid, nm_tid; // Note: assumes each tid has only one nm_tid
	float tscale_x, tscale_y, txoff, tyoff;
	bool emissive; // for lights

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale_x(1.0), tscale_y(1.0), txoff(0.0), tyoff(0.0), emissive(0) {}
	tid_nm_pair_t(int tid_, float txy=1.0) : tid(tid_), nm_tid(FLAT_NMAP_TEX), tscale_x(txy), tscale_y(txy), txoff(0.0), tyoff(0.0), emissive(0) {} // non-normal mapped 1:1 texture AR
	tid_nm_pair_t(int tid_, int nm_tid_, float tx, float ty, float xo=0.0, float yo=0.0) : tid(tid_), nm_tid(nm_tid_), tscale_x(tx), tscale_y(ty), txoff(xo), tyoff(yo), emissive(0) {}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}
	
	bool operator==(tid_nm_pair_t const &t) const {
		return (tid == t.tid && nm_tid == t.nm_tid && tscale_x == t.tscale_x && tscale_y == t.tscale_y && txoff == t.txoff && tyoff == t.tyoff && emissive == t.emissive);
	}
	bool operator!=(tid_nm_pair_t const &t) const {return !operator==(t);}
	int get_nm_tid() const {return ((nm_tid < 0) ? FLAT_NMAP_TEX : nm_tid);}
	colorRGBA get_avg_color() const {return texture_color(tid);}
	tid_nm_pair_t get_scaled_version(float scale) const {return tid_nm_pair_t(tid, nm_tid, scale*tscale_x, scale*tscale_y);}
	void set_gl(shader_t &s) const;
	void unset_gl(shader_t &s) const;
	void toggle_transparent_windows_mode();
};

struct building_tex_params_t {
	tid_nm_pair_t side_tex, roof_tex; // exterior
	tid_nm_pair_t wall_tex, ceil_tex, floor_tex, house_ceil_tex, house_floor_tex; // interior

	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0 || wall_tex.nm_tid >= 0 ||
		ceil_tex.nm_tid >= 0 || floor_tex.nm_tid >= 0 || house_ceil_tex.nm_tid >= 0 || house_floor_tex.nm_tid >= 0);}
};

struct color_range_t {

	float grayscale_rand;
	colorRGBA cmin, cmax; // alpha is unused?

	color_range_t() : grayscale_rand(0.0), cmin(WHITE), cmax(WHITE) {}
	void gen_color(colorRGBA &color, rand_gen_t &rgen) const;
};

struct building_mat_t : public building_tex_params_t {

	bool no_city, add_windows, add_wind_lights;
	unsigned min_levels, max_levels, min_sides, max_sides;
	float place_radius, max_delta_z, max_rot_angle, min_level_height, min_alt, max_alt, house_prob, house_scale_min, house_scale_max;
	float split_prob, cube_prob, round_prob, asf_prob, min_fsa, max_fsa, min_asf, max_asf, wind_xscale, wind_yscale, wind_xoff, wind_yoff;
	cube_t pos_range, prev_pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color; // exterior
	colorRGBA window_color, wall_color, ceil_color, floor_color, house_ceil_color, house_floor_color;

	building_mat_t() : no_city(0), add_windows(0), add_wind_lights(0), min_levels(1), max_levels(1), min_sides(4), max_sides(4), place_radius(0.0),
		max_delta_z(0.0), max_rot_angle(0.0), min_level_height(0.0), min_alt(-1000), max_alt(1000), house_prob(0.0), house_scale_min(1.0), house_scale_max(1.0),
		split_prob(0.0), cube_prob(1.0), round_prob(0.0), asf_prob(0.0), min_fsa(0.0), max_fsa(0.0), min_asf(0.0), max_asf(0.0), wind_xscale(1.0),
		wind_yscale(1.0), wind_xoff(0.0), wind_yoff(0.0), pos_range(-100,100,-100,100,0,0), prev_pos_range(all_zeros), sz_range(1,1,1,1,1,1),
		window_color(GRAY), wall_color(WHITE), ceil_color(WHITE), floor_color(LT_GRAY), house_ceil_color(WHITE), house_floor_color(WHITE) {}
	float gen_size_scale(rand_gen_t &rgen) const {return ((house_scale_min == house_scale_max) ? house_scale_min : rgen.rand_uniform(house_scale_min, house_scale_max));}
	void update_range(vector3d const &range_translate);
	void set_pos_range(cube_t const &new_pos_range) {prev_pos_range = pos_range; pos_range = new_pos_range;}
	void restore_prev_pos_range() {
		if (!prev_pos_range.is_all_zeros()) {pos_range = prev_pos_range;}
	}
	float get_window_tx() const;
	float get_window_ty() const;
	float get_floor_spacing() const {return 1.0/(2.0*get_window_ty());}
};

struct building_params_t {

	bool flatten_mesh, has_normal_map, tex_mirror, tex_inv_y, tt_only, infinite_buildings, dome_roof, onion_roof, enable_people_ai;
	unsigned num_place, num_tries, cur_prob;
	float ao_factor, sec_extra_spacing;
	float window_width, window_height, window_xspace, window_yspace; // windows
	float wall_split_thresh; // interiors
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix, mat_gen_ix_city, mat_gen_ix_nocity; // {any, city_only, non_city}
	vector<unsigned> rug_tids, picture_tids, sheet_tids;

	building_params_t(unsigned num=0) : flatten_mesh(0), has_normal_map(0), tex_mirror(0), tex_inv_y(0), tt_only(0), infinite_buildings(0), dome_roof(0),
		onion_roof(0), enable_people_ai(0), num_place(num), num_tries(10), cur_prob(1), ao_factor(0.0), sec_extra_spacing(0.0), window_width(0.0), window_height(0.0),
		window_xspace(0.0), window_yspace(0.0), wall_split_thresh(4.0), range_translate(zero_vector) {}
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
	bool windows_enabled  () const {return (window_width > 0.0 && window_height > 0.0 && window_xspace > 0.0 && window_yspace);} // all must be specified as nonzero
	bool gen_inf_buildings() const {return (infinite_buildings && world_mode == WMODE_INF_TERRAIN);}
	float get_window_width_fract () const {assert(windows_enabled()); return window_width /(window_width  + window_xspace);}
	float get_window_height_fract() const {assert(windows_enabled()); return window_height/(window_height + window_yspace);}
	float get_window_tx() const {assert(windows_enabled()); return 1.0f/(window_width  + window_xspace);}
	float get_window_ty() const {assert(windows_enabled()); return 1.0f/(window_height + window_yspace);}
	void add_cur_mat();

	void finalize() {
		if (materials.empty()) {add_cur_mat();} // add current (maybe default) material
	}
	building_mat_t const &get_material(unsigned mat_ix) const {
		assert(mat_ix < materials.size());
		return materials[mat_ix];
	}
	vector<unsigned> const &get_mat_list(bool city_only, bool non_city_only) const {
		return (city_only ? mat_gen_ix_city : (non_city_only ? mat_gen_ix_nocity : mat_gen_ix));
	}
	unsigned choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only) const;
	void set_pos_range(cube_t const &pos_range);
	void restore_prev_pos_range();
};

class building_draw_t;

struct building_geom_t { // describes the physical shape of a building
	unsigned num_sides;
	uint8_t door_sides[4]; // bit mask for 4 door sides, one per base part
	bool half_offset, is_pointed;
	float rot_sin, rot_cos, flat_side_amt, alt_step_factor, start_angle; // rotation in XY plane, around Z (up) axis
	//float roof_recess;

	building_geom_t(unsigned ns=4, float rs=0.0, float rc=1.0, bool pointed=0) : num_sides(ns), half_offset(0), is_pointed(pointed),
		rot_sin(rs), rot_cos(rc), flat_side_amt(0.0), alt_step_factor(0.0), start_angle(0.0)
	{
		door_sides[0] = door_sides[1] = door_sides[2] = door_sides[3] = 0;
	}
	bool is_rotated() const {return (rot_sin != 0.0);}
	bool is_cube()    const {return (num_sides == 4);}
	bool is_simple_cube()    const {return (is_cube() && !half_offset && flat_side_amt == 0.0 && alt_step_factor == 0.0);}
	bool use_cylinder_coll() const {return (num_sides > 8 && flat_side_amt == 0.0);} // use cylinder collision if not a cube, triangle, octagon, etc. (approximate)
};

struct tquad_with_ix_t : public tquad_t {
	// roof, roof access cover, wall, chimney cap, house door, building door, garage door, interior door, roof door
	enum {TYPE_ROOF=0, TYPE_ROOF_ACC, TYPE_WALL, TYPE_CCAP, TYPE_HDOOR, TYPE_BDOOR, TYPE_GDOOR, TYPE_IDOOR, TYPE_IDOOR2, TYPE_RDOOR};
	bool is_exterior_door() const {return (type == TYPE_HDOOR || type == TYPE_BDOOR || type == TYPE_GDOOR || type == TYPE_RDOOR);}
	bool is_interior_door() const {return (type == TYPE_IDOOR || type == TYPE_IDOOR2);}

	unsigned type;
	tquad_with_ix_t(unsigned npts_=0) : tquad_t(npts_), type(TYPE_ROOF) {}
	tquad_with_ix_t(tquad_t const &t, unsigned type_) : tquad_t(t), type(type_) {}
};

struct vertex_range_t {
	int draw_ix; // -1 is unset
	unsigned start, end;
	vertex_range_t() : draw_ix(-1), start(0), end(0) {}
	vertex_range_t(unsigned s, unsigned e, int ix=-1) : draw_ix(ix), start(s), end(e) {}
};

struct draw_range_t {
	// intended for building interiors, which don't have many materials; may need to increase MAX_DRAW_BLOCKS later
	vertex_range_t vr[MAX_DRAW_BLOCKS]; // quad verts only for now
};

enum room_object    {TYPE_NONE =0, TYPE_TABLE, TYPE_CHAIR, TYPE_STAIR, TYPE_ELEVATOR, TYPE_LIGHT, TYPE_RUG, TYPE_PICTURE, TYPE_WBOARD, TYPE_BOOK,
	                 TYPE_BCASE, TYPE_TCAN, TYPE_DESK, TYPE_BED, TYPE_WINDOW, TYPE_TOILET, TYPE_SINK, TYPE_TUB, TYPE_FRIDGE, TYPE_STOVE,
	                 TYPE_TV, TYPE_COUCH, NUM_TYPES};
enum room_obj_shape {SHAPE_CUBE=0, SHAPE_CYLIN, SHAPE_STAIRS_U, SHAPE_TALL};
enum room_type      {RTYPE_NOTSET=0, RTYPE_HALL, RTYPE_STAIRS, RTYPE_OFFICE, RTYPE_BATH, RTYPE_BED, RTYPE_KITCHEN, RTYPE_LIVING, RTYPE_DINING, RTYPE_STUDY, RTYPE_GARAGE, RTYPE_SHED, NUM_RTYPES};
std::string const room_names[NUM_RTYPES] = {"Not Set", "Hallway", "Stairs", "Office", "Bathroom", "Bedroom", "Kitchen", "Living Room", "Dining Room", "Study", "Garage", "Shed"};
enum stairs_shape   {SHAPE_STRAIGHT=0, SHAPE_U, SHAPE_WALLED};
enum {ROOM_WALL_INT=0, ROOM_WALL_SEP, ROOM_WALL_EXT};
enum {OBJ_MODEL_TOILET=0, OBJ_MODEL_SINK, OBJ_MODEL_TUB, OBJ_MODEL_FRIDGE, OBJ_MODEL_STOVE, OBJ_MODEL_TV, OBJ_MODEL_COUCH, NUM_OBJ_MODELS};

// object flags, currently used for room lights
uint8_t const RO_FLAG_LIT     = 0x01; // light is on
uint8_t const RO_FLAG_TOS     = 0x02; // at top of stairs
uint8_t const RO_FLAG_RSTAIRS = 0x04; // in a room with stairs
uint8_t const RO_FLAG_INVIS   = 0x08; // invisible
uint8_t const RO_FLAG_NOCOLL  = 0x10; // no collision detection
uint8_t const RO_FLAG_OPEN    = 0x20; // open, for elevators and maybe eventually doors
uint8_t const RO_FLAG_NODYNAM = 0x40; // for light shadow maps

struct room_object_t : public cube_t {
	bool dim, dir;
	uint8_t flags, room_id; // for at most 256 rooms per floor
	uint16_t obj_id; // currently only used for lights and random property hashing
	room_object type;
	room_obj_shape shape;
	float light_amt;
	colorRGBA color;

	room_object_t() : dim(0), dir(0), flags(0), room_id(0), obj_id(0), type(TYPE_NONE), shape(SHAPE_CUBE), light_amt(1.0) {}
	room_object_t(cube_t const &c, room_object type_, uint8_t rid, bool dim_=0, bool dir_=0, uint8_t f=0, float light=1.0,
		room_obj_shape shape_=room_obj_shape::SHAPE_CUBE, colorRGBA const color_=WHITE) :
		cube_t(c), dim(dim_), dir(dir_), flags(f), room_id(rid), obj_id(0), type(type_), shape(shape_), light_amt(light), color(color_)
	{assert(is_strictly_normalized());}
	bool is_lit    () const {return (flags & RO_FLAG_LIT);}
	bool has_stairs() const {return (flags & (RO_FLAG_TOS | RO_FLAG_RSTAIRS));}
	bool is_visible() const {return !(flags & RO_FLAG_INVIS);}
	bool no_coll   () const {return (flags & RO_FLAG_NOCOLL);}
	void toggle_lit_state() {flags ^= RO_FLAG_LIT;}
	static bool enable_rugs();
	static bool enable_pictures();
	int get_rug_tid() const;
	int get_picture_tid() const;
	int get_sheet_tid() const;
	colorRGBA get_color() const;
};

struct rgeom_storage_t {
	typedef vert_norm_comp_tc_color vertex_t;
	vector<vertex_t> quad_verts, itri_verts;
	vector<unsigned> indices; // for indexed_quad_verts
	tid_nm_pair_t tex; // used sort of like a map key

	rgeom_storage_t() {}
	rgeom_storage_t(tid_nm_pair_t const &tex_) : tex(tex_) {}
	void clear();
	void swap_vectors(rgeom_storage_t &s);
	void swap(rgeom_storage_t &s);
	unsigned get_tot_vert_capacity() const {return (quad_verts.capacity() + itri_verts.capacity());}
};

class rgeom_mat_t : public rgeom_storage_t { // simplified version of building_draw_t::draw_block_t

	vbo_wrap_t vbo;
	unsigned ivbo;
public:
	unsigned num_qverts, num_itverts, num_ixs; // for drawing
	bool en_shadows;

	rgeom_mat_t(tid_nm_pair_t const &tex_) : rgeom_storage_t(tex_), ivbo(0), num_qverts(0), num_itverts(0), num_ixs(0), en_shadows(0) {}
	unsigned get_tot_vert_count() const {return (num_qverts + num_itverts);}
	void enable_shadows() {en_shadows = 1;}
	void clear();
	void add_cube_to_verts(cube_t const &c, colorRGBA const &color, vector3d const &tex_origin, unsigned skip_faces=0, bool swap_tex_st=0, bool mirror_x=0, bool mirror_y=0);
	void add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided=0, bool ts_tb=0, bool inv_tb=0, float rs_bot=1.0, float rs_top=1.0);
	void create_vbo();
	void draw(shader_t &s, bool shadow_only);
};

struct building_materials_t : public vector<rgeom_mat_t> {
	void clear();
	unsigned count_all_verts() const;
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows);
	void create_vbos();
	void draw(shader_t &s, bool shadow_only);
};

struct obj_model_inst_t {
	unsigned obj_id, model_id;
	colorRGBA color;
	obj_model_inst_t(unsigned oid, unsigned mid, colorRGBA const &c) : obj_id(oid), model_id(mid), color(c) {}
};

struct building_room_geom_t {

	bool has_elevators, has_pictures;
	unsigned char num_pic_tids;
	float obj_scale;
	unsigned stairs_start; // index of first object of TYPE_STAIR
	vector3d tex_origin;
	vector<room_object_t> objs; // for drawing and collision detection
	vector<obj_model_inst_t> obj_model_insts;
	building_materials_t mats_static, mats_small, mats_dynamic; // {large static, small static, dynamic} materials
	vect_cube_t light_bcubes;

	building_room_geom_t(vector3d const &tex_origin_) : has_elevators(0), has_pictures(0), num_pic_tids(0), obj_scale(1.0), stairs_start(0), tex_origin(tex_origin_) {}
	bool empty() const {return objs.empty();}
	void clear();
	void clear_materials();
	unsigned get_num_verts() const {return (mats_static.count_all_verts() + mats_small.count_all_verts() + mats_dynamic.count_all_verts());}
	rgeom_mat_t &get_material(tid_nm_pair_t const &tex, bool inc_shadows=0, bool dynamic=0, bool small=0);
	rgeom_mat_t &get_wood_material(float tscale);
	void add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale);
	void add_table(room_object_t const &c, float tscale);
	void add_chair(room_object_t const &c, float tscale);
	void add_stair(room_object_t const &c, float tscale, vector3d const &tex_origin);
	void add_elevator(room_object_t const &c, float tscale);
	void add_light(room_object_t const &c, float tscale);
	void add_rug(room_object_t const &c);
	void add_picture(room_object_t const &c);
	void add_book_title(std::string const &title, cube_t const &title_area, rgeom_mat_t &mat, colorRGBA const &color, unsigned hdim, unsigned tdim, unsigned wdim, bool cdir, bool ldir, bool wdir);
	void add_book(room_object_t const &c, bool inc_lg, bool inc_sm, float tilt_angle=0.0, unsigned extra_skip_faces=0, bool no_title=0);
	void add_bookcase(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale, bool no_shelves=0, float sides_scale=1.0);
	void add_desk(room_object_t const &c, float tscale);
	void add_bed(room_object_t const &c, bool inc_lg, bool inc_sm, float tscale);
	void add_window(room_object_t const &c, float tscale);
	void add_tub_outer(room_object_t const &c);
	void add_tv_picture(room_object_t const &c);
	void add_trashcan(room_object_t const &c);
	void create_static_vbos(bool small_objs);
	void create_dynamic_vbos();
	void draw(shader_t &s, vector3d const &xlate, bool shadow_only, bool inc_small);
};

struct elevator_t : public cube_t {
	bool dim, dir, open, at_edge; // door dim/dir
	elevator_t() : dim(0), dir(0), open(0), at_edge(0) {}
	elevator_t(cube_t const &c, bool dim_, bool dir_, bool open_, bool at_edge_) : cube_t(c), dim(dim_), dir(dir_), open(open_), at_edge(at_edge_) {assert(is_strictly_normalized());}
	float get_wall_thickness() const {return 0.02*get_sz_dim(!dim);}
	float get_frame_width   () const {return 0.20*get_sz_dim(!dim);}
	unsigned get_door_face_id() const {return (2*dim + dir);}
	unsigned get_coll_cubes(cube_t cubes[5]) const; // returns 1 or 5 cubes
};

struct room_t : public cube_t {
	bool has_stairs, has_elevator, no_geom, is_hallway, is_office, is_sec_bldg; // or room_type index?
	uint8_t ext_sides; // sides that have exteriors, and likely windows (bits for x1, x2, y1, y2)
	//uint8_t sides_with_doors; // is this useful/needed?
	uint8_t part_id, num_lights;
	room_type rtype; // this applies to the first floor because some rooms can have variable per-floor assignment
	uint64_t lit_by_floor;
	room_t() : has_stairs(0), has_elevator(0), no_geom(0), is_hallway(0), is_office(0), is_sec_bldg(0), ext_sides(0), part_id(0), num_lights(0), rtype(RTYPE_NOTSET), lit_by_floor(0) {}
	room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_);
	void assign_to(room_type rt, unsigned floor=0);
	float get_light_amt() const;
};

struct landing_t : public cube_t {
	bool for_elevator, dim, dir;
	uint8_t floor;
	stairs_shape shape;
	landing_t() : for_elevator(0), dim(0), dir(0), floor(0), shape(SHAPE_STRAIGHT) {}
	landing_t(cube_t const &c, bool e, uint8_t f, bool dim_, bool dir_, stairs_shape shape_=SHAPE_STRAIGHT) :
		cube_t(c), for_elevator(e), dim(dim_), dir(dir_), floor(f), shape(shape_) {assert(is_strictly_normalized());}
	unsigned get_face_id() const {return (2*dim + dir);}
};

struct door_t : public cube_t {
	bool dim, open_dir, open;
	door_t() : dim(0), open_dir(0), open(0) {}
	door_t(cube_t const &c, bool dim_, bool dir, bool open_) : cube_t(c), dim(dim_), open_dir(dir), open(open_) {assert(is_strictly_normalized());}
};
typedef vector<door_t> vect_door_t;

enum {ROOF_OBJ_BLOCK=0, ROOF_OBJ_ANT, ROOF_OBJ_WALL, ROOF_OBJ_ECAP, ROOF_OBJ_AC, ROOF_OBJ_SCAP};
enum {ROOF_TYPE_FLAT=0, ROOF_TYPE_SLOPE, ROOF_TYPE_PEAK, ROOF_TYPE_DOME, ROOF_TYPE_ONION};

struct roof_obj_t : public cube_t {
	uint8_t type;
	roof_obj_t(uint8_t type_=ROOF_OBJ_BLOCK) : type(type_) {}
	roof_obj_t(cube_t const &c, uint8_t type_=ROOF_OBJ_BLOCK) : cube_t(c), type(type_) {assert(is_strictly_normalized());}
};
typedef vector<roof_obj_t> vect_roof_obj_t;

struct stairwell_t : public cube_t {
	stairs_shape shape;
	uint8_t num_floors;
	bool dim, dir, roof_access, stack_conn;
	stairwell_t() : shape(SHAPE_STRAIGHT), num_floors(0), dim(0), dir(0), roof_access(0), stack_conn(0) {}
	stairwell_t(cube_t const &c, unsigned n, bool dim_, bool dir_, stairs_shape s=SHAPE_STRAIGHT, bool r=0, bool sc=0) :
		cube_t(c), shape(s), num_floors(n), dim(dim_), dir(dir_), roof_access(r), stack_conn(sc) {}
};
typedef vector<stairwell_t> vect_stairwell_t;

// may as well make this its own class, since it could get large and it won't be used for every building
struct building_interior_t {
	vect_cube_t floors, ceilings, walls[2]; // walls are split by dim
	vect_stairwell_t stairwells;
	vector<door_t> doors;
	vector<landing_t> landings; // for stairs and elevators
	vector<room_t> rooms;
	vector<elevator_t> elevators;
	std::unique_ptr<building_room_geom_t> room_geom;
	std::unique_ptr<building_nav_graph_t> nav_graph;
	draw_range_t draw_range;
	uint64_t top_ceilings_mask; // bit mask for ceilings that are on the top floor and have no floor above them

	building_interior_t();
	~building_interior_t();
	bool is_cube_close_to_doorway(cube_t const &c, float dmin=0.0f, bool inc_open=0) const;
	bool is_blocked_by_stairs_or_elevator(cube_t const &c, float dmin=0.0f, bool elevators_only=0) const;
	bool is_blocked_by_stairs_or_elevator_no_expand(cube_t const &c, float dmin=0.0f) const;
	void finalize();
	bool update_elevators(point const &player_pos, float floor_thickness);
	void get_avoid_cubes(vect_cube_t &avoid, float z1, float z2) const;
};

struct building_stats_t {
	unsigned nbuildings, nparts, ndetails, ntquads, ndoors, ninterior, nrooms, nceils, nfloors, nwalls, nrgeom, nobjs, nverts;
	building_stats_t() : nbuildings(0), nparts(0), ndetails(0), ntquads(0), ndoors(0), ninterior(0), nrooms(0), nceils(0), nfloors(0), nwalls(0), nrgeom(0), nobjs(0), nverts(0) {}
};

struct building_loc_t {
	int part_ix, room_ix, stairs_ix; // -1 is not contained; what about elevator_ix?
	unsigned floor;
	building_loc_t() : part_ix(-1), room_ix(-1), stairs_ix(-1), floor(0) {}
	bool operator==(building_loc_t const &loc) const {return (part_ix == loc.part_ix && room_ix == loc.room_ix && stairs_ix == loc.stairs_ix && floor == loc.floor);}
};

enum {AI_STOP=0, AI_WAITING, AI_NEXT_PT, AI_BEGIN_PATH, AI_AT_DEST, AI_MOVING};

struct building_ai_state_t {
	bool is_first_path;
	unsigned cur_room, dest_room; // Note: cur_room and dest_room may not be needed
	vector<point> path; // stored backwards, next point on path is path.back()

	building_ai_state_t() : is_first_path(1), cur_room(0), dest_room(0) {}
	void next_path_pt(pedestrian_t &person, bool same_floor);
};

struct colored_cube_t;
typedef vector<colored_cube_t> vect_colored_cube_t;
class cube_bvh_t;
class building_indir_light_mgr_t;

struct building_t : public building_geom_t {

	unsigned mat_ix;
	uint8_t hallway_dim, real_num_parts, roof_type; // main hallway dim: 0=x, 1=y, 2=none
	bool is_house, has_chimney, has_garage, has_shed, has_courtyard, has_complex_floorplan;
	colorRGBA side_color, roof_color, detail_color;
	cube_t bcube, pri_hall;
	vect_cube_t parts;
	vect_roof_obj_t details; // cubes on the roof - antennas, AC units, etc.
	vector<tquad_with_ix_t> roof_tquads, doors;
	std::shared_ptr<building_interior_t> interior;
	vertex_range_t ext_side_qv_range;
	point tree_pos; // (0,0,0) is unplaced/no tree
	float ao_bcz2;

	friend class building_indir_light_mgr_t;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_), hallway_dim(2), real_num_parts(0), roof_type(ROOF_TYPE_FLAT), is_house(0), has_chimney(0),
		has_garage(0), has_shed(0), has_courtyard(0), has_complex_floorplan(0), side_color(WHITE), roof_color(WHITE), detail_color(BLACK), ao_bcz2(0.0) {}
	building_t(building_geom_t const &bg) : building_geom_t(bg), mat_ix(0), hallway_dim(2), real_num_parts(0), roof_type(ROOF_TYPE_FLAT),
		is_house(0), has_chimney(0), has_garage(0), has_shed(0), has_courtyard(0), has_complex_floorplan(0), ao_bcz2(0.0) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	bool has_interior () const {return bool(interior);}
	bool has_room_geom() const {return (has_interior() && interior->room_geom);}
	bool has_sec_bldg () const {return (has_garage || has_shed);}
	bool has_pri_hall () const {return (hallway_dim <= 1);} // otherswise == 2
	colorRGBA get_avg_side_color  () const {return side_color  .modulate_with(get_material().side_tex.get_avg_color());}
	colorRGBA get_avg_roof_color  () const {return roof_color  .modulate_with(get_material().roof_tex.get_avg_color());}
	colorRGBA get_avg_detail_color() const {return detail_color.modulate_with(get_material().roof_tex.get_avg_color());}
	building_mat_t const &get_material() const;
	float get_window_vspace  () const {return get_material().get_floor_spacing();}
	float get_floor_thickness() const {return FLOOR_THICK_VAL*get_window_vspace();}
	float get_wall_thickness () const {return 0.5*get_floor_thickness();}
	float get_door_height    () const {return 0.9f*get_window_vspace();} // set height based on window spacing, 90% of a floor height (may be too large)
	unsigned get_real_num_parts() const {return (is_house ? min(2U, unsigned(parts.size() - has_chimney)) : parts.size());}
	void gen_rotation(rand_gen_t &rgen);
	void set_z_range(float z1, float z2);
	bool check_part_contains_pt_xy(cube_t const &part, point const &pt, vector<point> &points) const;
	bool check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const;
	vect_cube_t::const_iterator get_real_parts_end() const {return (parts.begin() + real_num_parts);}
	vect_cube_t::const_iterator get_real_parts_end_inc_sec() const {return (get_real_parts_end() + has_sec_bldg());}
	cube_t const &get_sec_bldg() const {assert(has_sec_bldg()); assert(real_num_parts < parts.size()); return parts[real_num_parts];}
	void end_add_parts() {assert(parts.size() < 256); real_num_parts = uint8_t(parts.size());}

	bool check_sphere_coll(point const &pos, float radius, bool xy_only, vector<point> &points, vector3d *cnorm=nullptr) const {
		point pos2(pos);
		return check_sphere_coll(pos2, pos, vect_cube_t(), zero_vector, radius, xy_only, points, cnorm);
	}
	bool check_sphere_coll(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, vector3d const &xlate, float radius, bool xy_only,
		vector<point> &points, vector3d *cnorm=nullptr, bool check_interior=0) const;
	bool check_sphere_coll_interior(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, float radius, bool xy_only, vector3d *cnorm=nullptr) const;
	unsigned check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t, vector<point> &points, bool occlusion_only=0, bool ret_any_pt=0, bool no_coll_pt=0) const;
	bool check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const;
	bool ray_cast_interior(point const &pos, vector3d const &dir, cube_bvh_t const &bvh, point &cpos, vector3d &cnorm, colorRGBA &ccolor) const;
	void create_building_volume_light_texture(unsigned bix, point const &target, unsigned &tid) const;
	bool ray_cast_camera_dir(point const &camera_bs, point &cpos, colorRGBA &ccolor) const;
	void calc_bcube_from_parts();
	void adjust_part_zvals_for_floor_spacing(cube_t &c) const;
	void gen_geometry(int rseed1, int rseed2);
	cube_t place_door(cube_t const &base, bool dim, bool dir, float door_height, float door_center, float door_pos, float door_center_shift, float width_scale, bool can_fail, rand_gen_t &rgen);
	void gen_house(cube_t const &base, rand_gen_t &rgen);
	bool add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_building, bool roof_access=0);
	float gen_peaked_roof(cube_t const &top_, float peak_height, bool dim, float extend_to, float max_dz, unsigned skip_side_tri);
	float gen_hipped_roof(cube_t const &top_, float peak_height, float extend_to);
	void place_roof_ac_units(unsigned num, float sz_scale, cube_t const &bounds, vect_cube_t const &avoid, bool avoid_center, rand_gen_t &rgen);
	void add_roof_walls(cube_t const &c, float wall_width, bool overlap_corners, cube_t out[4]);
	void gen_details(rand_gen_t &rgen, bool is_rectangle);
	int get_num_windows_on_side(float xy1, float xy2) const;
	bool interior_enabled() const;
	void gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes);
	void add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned part_ix, unsigned num_floors,
		unsigned rooms_start, bool use_hallway, bool first_part_this_stack, float window_hspacing[2], float window_border);
	void connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part);
	void gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes);
	void add_stairs_and_elevators(rand_gen_t &rgen);
	void gen_building_doors_if_needed(rand_gen_t &rgen);
	void maybe_add_special_roof(rand_gen_t &rgen);
	void gen_sloped_roof(rand_gen_t &rgen, cube_t const &top);
	void add_roof_to_bcube();
	void gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax);
	void get_all_drawn_verts(building_draw_t &bdraw, bool get_exterior, bool get_interior);
	void get_all_drawn_window_verts(building_draw_t &bdraw, bool lights_pass, float offset_scale=1.0, point const *const only_cont_pt=nullptr) const;
	bool get_nearby_ext_door_verts(building_draw_t &bdraw, shader_t &s, point const &pos, float dist, unsigned &door_type) const;
	bool find_door_close_to_point(tquad_with_ix_t &door, point const &pos, float dist) const;
	void get_split_int_window_wall_verts(building_draw_t &bdraw_front, building_draw_t &bdraw_back, point const &only_cont_pt, bool make_all_front=0) const;
	void get_ext_wall_verts_no_sec(building_draw_t &bdraw) const;
	void add_room_lights(vector3d const &xlate, unsigned building_id, bool camera_in_building, int ped_ix, vect_cube_t &ped_bcubes, cube_t &lights_bcube);
	bool toggle_room_light(point const &closest_to);
	void draw_room_geom(shader_t &s, vector3d const &xlate, bool shadow_only, bool inc_small);
	void gen_and_draw_room_geom(shader_t &s, vector3d const &xlate, vect_cube_t &ped_bcubes, unsigned building_ix, int ped_ix, bool shadow_only, bool inc_small);
	void add_split_roof_shadow_quads(building_draw_t &bdraw) const;
	void clear_room_geom();
	bool place_person(point &ppos, float radius, rand_gen_t &rgen) const;
	void update_grass_exclude_at_pos(point const &pos, vector3d const &xlate) const;
	void update_stats(building_stats_t &s) const;
	void build_nav_graph() const;
	unsigned count_connected_room_components() const;
	bool is_room_adjacent_to_ext_door(cube_t const &room) const;
	point get_center_of_room(unsigned room_ix) const;
	bool choose_dest_room(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const;
	bool find_route_to_point(point const &from, point const &to, float radius, bool is_first_path, vector<point> &path) const;
	void find_nearest_stairs(point const &p1, point const &p2, vector<unsigned> &nearest_stairs, bool straight_only, int part_ix=-1) const;
	int ai_room_update(building_ai_state_t &state, rand_gen_t &rgen, vector<pedestrian_t> &people, float delta_dir, unsigned person_ix, bool stay_on_one_floor=1) const;
	void move_person_to_not_collide(pedestrian_t &person, pedestrian_t const &other, point const &new_pos, float rsum, float coll_dist) const;
	building_loc_t get_building_loc_for_pt(point const &pt) const;
	bool place_obj_along_wall(room_object type, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
		bool is_lit, cube_t const &place_area, unsigned objs_start, float front_clearance=0.0, unsigned pref_orient=4, bool pref_centered=0, colorRGBA const &color=WHITE);
	bool place_model_along_wall(unsigned model_id, room_object type, float height, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
		bool is_lit, cube_t const &place_area, unsigned objs_start, float front_clearance=0.0, unsigned pref_orient=4, bool pref_centered=0, colorRGBA const &color=WHITE);
	void update_elevators(point const &player_pos);
private:
	cube_t get_walkable_room_bounds(room_t const &room) const;
	void get_exclude_cube(point const &pos, cube_t const &skip, cube_t &exclude) const;
	void add_door_to_bdraw(cube_t const &D, building_draw_t &bdraw, uint8_t door_type, bool dim, bool dir, bool opened, bool opens_out, bool exterior) const;
	void move_door_to_other_side_of_wall(tquad_with_ix_t &door, float dist_mult, bool invert_normal) const;
	tquad_with_ix_t set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj,
		bool exterior, bool opened, bool opens_out, bool opens_up, bool swap_sides) const;
	void clip_door_to_interior(tquad_with_ix_t &door, bool clip_to_floor) const;
	cube_t get_part_containing_pt(point const &pt) const;
	bool is_cube_close_to_doorway(cube_t const &c, float dmin=0.0, bool inc_open=0) const;
	bool is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, float room_pad=0.0) const;
	bool check_cube_intersect_walls(cube_t const &c) const;
	bool check_cube_contained_in_part(cube_t const &c) const;
	bool is_valid_stairs_elevator_placement(cube_t const &c, float door_pad, float stairs_pad, bool check_walls=1) const;
	bool clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const;
	void add_room(cube_t const &room, unsigned part_id, unsigned num_lights, bool is_hallway, bool is_office, bool is_sec_bldg=0);
	void add_or_extend_elevator(elevator_t const &elevator, bool add);
	void remove_intersecting_roof_cubes(cube_t const &c);
	bool overlaps_other_room_obj(cube_t const &c, unsigned objs_start) const;
	int classify_room_wall(room_t const &room, float zval, bool dim, bool dir, bool ret_sep_if_part_int_part_ext) const;
	bool add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos,
		colorRGBA const &chair_color, bool dim, bool dir, float tot_light_amt, bool is_lit);
	bool add_table_and_chairs(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
		point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, bool is_lit);
	void add_trashcan_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start, bool check_last_obj);
	bool add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start);
	bool add_desk_to_room    (rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, colorRGBA const &chair_color, float zval, unsigned room_id, float tot_light_amt, bool is_lit);
	bool add_bed_to_room     (rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, float tot_light_amt, bool is_lit);
	bool add_bathroom_objs   (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start);
	bool add_kitchen_objs    (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start, bool allow_adj_ext_door);
	bool add_livingroom_objs (rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start);
	void place_book_on_obj   (rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, bool is_lit, bool use_dim_dir);
	void add_rug_to_room     (rand_gen_t &rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit);
	bool hang_pictures_in_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start);
	void add_bathroom_windows(room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit);
	bool can_be_bedroom_or_bathroom(room_t const &room, bool on_first_floor) const;
	bool can_be_bathroom(room_t const &room) const;
	unsigned count_num_int_doors(room_t const &room) const;
	bool check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const;
	void split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen);
	bool test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const;
	void gather_interior_cubes(vect_colored_cube_t &cc) const;
	void order_lights_by_priority(point const &target, vector<unsigned> &light_ids) const;
	bool is_light_occluded(point const &lpos, point const &camera_bs) const;
	void clip_ray_to_walls(point const &p1, point &p2) const;
	void refine_light_bcube(point const &lpos, float light_radius, cube_t &light_bcube) const;
};

struct vect_building_t : public vector<building_t> {
	void ai_room_update(vector<building_ai_state_t> &ai_state, vector<pedestrian_t> &people, float delta_dir, rand_gen_t &rgen) const;
};

struct building_draw_utils {
	static void calc_normals(building_geom_t const &bg, vector<vector3d> &nv, unsigned ndiv);
	static void calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, vector<point> &pts, float expand=0.0);
};

class city_lights_manager_t {
protected:
	cube_t lights_bcube;
	float light_radius_scale, dlight_add_thresh;
	bool prev_had_lights;
public:
	city_lights_manager_t() : lights_bcube(all_zeros), light_radius_scale(1.0), dlight_add_thresh(0.0), prev_had_lights(0) {}
	virtual ~city_lights_manager_t() {}
	cube_t get_lights_bcube() const {return lights_bcube;}
	void add_player_flashlight(float radius_scale);
	void tighten_light_bcube_bounds(vector<light_source> const &lights);
	void clamp_to_max_lights(vector3d const &xlate, vector<light_source> &lights);
	bool begin_lights_setup(vector3d const &xlate, float light_radius, vector<light_source> &lights);
	void finalize_lights(vector<light_source> &lights);
	void setup_shadow_maps(vector<light_source> &light_sources, point const &cpos);
	virtual bool enable_lights() const = 0;
};

struct building_place_t {
	point p;
	unsigned bix;
	building_place_t() : bix(0) {}
	building_place_t(point const &p_, unsigned bix_) : p(p_), bix(bix_) {}
	bool operator<(building_place_t const &p) const {return (bix < p.bix);} // only compare building index
};
typedef vector<building_place_t> vect_building_place_t;

inline void clip_low_high(float &t0, float &t1) {
	if (fabs(t0 - t1) < 0.5) {t0 = t1 = 0.0;} // too small to have a window
	else {t0 = round_fp(t0); t1 = round_fp(t1);} // Note: round() is much faster than nearbyint(), and round_fp() is faster than round()
}

template<typename T> bool has_bcube_int(cube_t const &bcube, vector<T> const &bcubes, bool inc_adj=1) { // T must derive from cube_t
	for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {
		if (inc_adj ? c->intersects(bcube) : c->intersects_no_adj(bcube)) return 1;
	}
	return 0;
}

void do_xy_rotate(float rot_sin, float rot_cos, point const &center, point &pos);
void do_xy_rotate_normal(float rot_sin, float rot_cos, point &n);
void get_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state);
cube_t get_building_lights_bcube();
template<typename T> bool has_bcube_int_xy(cube_t const &bcube, vector<T> const &bcubes, float pad_dist=0.0);
void add_building_interior_lights(point const &xlate, cube_t &lights_bcube);
unsigned calc_num_floors(cube_t const &c, float window_vspacing, float floor_thickness);
void set_wall_width(cube_t &wall, float pos, float half_thick, bool dim);
void subtract_cube_from_cube(cube_t const &c, cube_t const &s, vect_cube_t &out);
void subtract_cube_from_cube_inplace(cube_t const &s, vect_cube_t &cubes, unsigned &ix, unsigned &iter_end);
template<typename T> void subtract_cubes_from_cube(cube_t const &c, T const &sub, vect_cube_t &out, vect_cube_t &out2);
bool subtract_cube_from_cubes(cube_t const &s, vect_cube_t &cubes, vect_cube_t *holes=nullptr, bool clip_in_z=0);
int get_rect_panel_tid();
int get_bath_wind_tid ();
int get_normal_map_for_bldg_tid(int tid);
// functions in city_gen.cc
void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, int use_smap, int use_bmap,
	float min_alpha=0.0, bool force_tsl=0, float pcf_scale=1.0, bool use_texgen=0, bool indir_lighting=0);
void enable_animations_for_shader(shader_t &s);
void setup_city_lights(vector3d const &xlate);
void draw_peds_in_building(int first_ped_ix, unsigned bix, shader_t &s, vector3d const &xlate, bool dlight_shadow_only); // from city_gen.cpp
void get_ped_bcubes_for_building(int first_ped_ix, unsigned bix, vect_cube_t &bcubes); // from city_gen.cpp
vector3d get_nom_car_size();
void draw_cars_in_garages(vector3d const &xlate, bool shadow_only);

