// 3D World - City Generation
// by Frank Gennari
// 2/10/18

#include "3DWorld.h"
#include "mesh.h"
#include "heightmap.h"
#include "file_utils.h"
#include "draw_utils.h"
#include "shaders.h"
#include "model3d.h"
#include "lightmap.h"
#include "buildings.h"

using std::string;

bool const CHECK_HEIGHT_BORDER_ONLY = 1; // choose building site to minimize edge discontinuity rather than amount of land that needs to be modified
float const ROAD_HEIGHT             = 0.002;
float const OUTSIDE_TERRAIN_HEIGHT  = 0.0;
float const CAR_LANE_OFFSET         = 0.15; // in units of road width
float const CONN_ROAD_SPEED_MULT    = 2.0; // twice the speed limit on connector roads
float const PARK_SPACE_WIDTH        = 1.6;
float const PARK_SPACE_LENGTH       = 1.8;
float const STREETLIGHT_BEAMWIDTH   = 0.25;
float const CITY_LIGHT_FALLOFF      = 0.2;
vector3d const CAR_SIZE(0.30, 0.13, 0.08); // {length, width, height} in units of road width
unsigned  const CONN_CITY_IX((1<<16)-1); // uint16_t max

enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   TID_PARK_LOT,  NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, TYPE_PARK_LOT, NUM_RD_TYPES};
enum {TURN_NONE=0, TURN_LEFT, TURN_RIGHT, TURN_UNSPEC};
unsigned const CONN_TYPE_NONE = 0;
colorRGBA const stoplight_colors[3] = {GREEN, YELLOW, RED};
colorRGBA const road_colors[NUM_RD_TYPES] = {WHITE, WHITE, WHITE, WHITE, WHITE, WHITE}; // parking lots are darker than roads


extern bool enable_dlight_shadows, dl_smap_enabled;
extern int rand_gen_index, display_mode, animate2;
extern unsigned shadow_map_sz;
extern float water_plane_z, shadow_map_pcf_offset, cobj_z_bias, fticks;
extern double tfticks;
extern vector<light_source> dl_sources;


void add_dynamic_lights_city(cube_t const &scene_bcube);


struct city_params_t {

	unsigned num_cities, num_samples, num_conn_tries, city_size_min, city_size_max, city_border, road_border, slope_width;
	float road_width, road_spacing, conn_road_seg_len, max_road_slope;
	// cars
	unsigned num_cars;
	float car_speed;
	// parking lots
	unsigned min_park_spaces, min_park_rows;
	float min_park_density, max_park_density;
	// lighting
	bool car_shadows;
	unsigned max_lights, max_shadow_maps;

	city_params_t() : num_cities(0), num_samples(100), num_conn_tries(50), city_size_min(0), city_size_max(0), city_border(0), road_border(0),
		slope_width(0), road_width(0.0), road_spacing(0.0), conn_road_seg_len(1000.0), max_road_slope(1.0), num_cars(0), car_speed(0.0),
		min_park_spaces(12), min_park_rows(1), min_park_density(0.0), max_park_density(1.0), car_shadows(0), max_lights(1024), max_shadow_maps(0) {}
	bool enabled() const {return (num_cities > 0 && city_size_min > 0);}
	bool roads_enabled() const {return (road_width > 0.0 && road_spacing > 0.0);}
	float get_road_ar() const {return nearbyint(road_spacing/road_width);} // round to nearest texture multiple
	static bool read_error(string const &str) {cout << "Error reading city config option " << str << "." << endl; return 0;}

	bool read_option(FILE *fp) {
		char strc[MAX_CHARS] = {0};
		if (!read_str(fp, strc)) return 0;
		string const str(strc);

		if (str == "num_cities") {
			if (!read_uint(fp, num_cities)) {return read_error(str);}
		}
		else if (str == "num_samples") {
			if (!read_uint(fp, num_samples) || num_samples == 0) {return read_error(str);}
		}
		else if (str == "num_conn_tries") {
			if (!read_uint(fp, num_conn_tries) || num_conn_tries == 0) {return read_error(str);}
		}
		else if (str == "city_size_min") {
			if (!read_uint(fp, city_size_min)) {return read_error(str);}
			if (city_size_max == 0) {city_size_max = city_size_min;}
			if (city_size_max < city_size_min) {return read_error(str);}
		}
		else if (str == "city_size_max") {
			if (!read_uint(fp, city_size_max)) {return read_error(str);}
			if (city_size_min == 0) {city_size_min = city_size_max;}
			if (city_size_max < city_size_min) {return read_error(str);}
		}
		else if (str == "city_border") {
			if (!read_uint(fp, city_border)) {return read_error(str);}
		}
		else if (str == "road_border") {
			if (!read_uint(fp, road_border)) {return read_error(str);}
		}
		else if (str == "slope_width") {
			if (!read_uint(fp, slope_width)) {return read_error(str);}
		}
		else if (str == "road_width") {
			if (!read_float(fp, road_width) || road_width < 0.0) {return read_error(str);}
		}
		else if (str == "road_spacing") {
			if (!read_float(fp, road_spacing) || road_spacing < 0.0) {return read_error(str);}
		}
		else if (str == "conn_road_seg_len") {
			if (!read_float(fp, conn_road_seg_len) || conn_road_seg_len <= 0.0) {return read_error(str);}
		}
		else if (str == "max_road_slope") {
			if (!read_float(fp, max_road_slope) || max_road_slope <= 0.0) {return read_error(str);}
		}
		// cars
		else if (str == "num_cars") {
			if (!read_uint(fp, num_cars)) {return read_error(str);}
		}
		else if (str == "car_speed") {
			if (!read_float(fp, car_speed) || car_speed < 0.0) {return read_error(str);}
		}
		// parking lots
		else if (str == "min_park_spaces") { // with default road parameters, can be up to 28
			if (!read_uint(fp, min_park_spaces)) {return read_error(str);}
		}
		else if (str == "min_park_rows") { // with default road parameters, can be up to 8
			if (!read_uint(fp, min_park_rows)) {return read_error(str);}
		}
		else if (str == "min_park_density") {
			if (!read_float(fp, min_park_density)) {return read_error(str);}
		}
		else if (str == "max_park_density") {
			if (!read_float(fp, max_park_density) || max_park_density < 0.0) {return read_error(str);}
		}
		// lighting
		else if (str == "max_lights") {
			if (!read_uint(fp, max_lights)) {return read_error(str);}
		}
		else if (str == "max_shadow_maps") {
			if (!read_uint(fp, max_shadow_maps)) {return read_error(str);}
		}
		else if (str == "car_shadows") {
			if (!read_bool(fp, car_shadows)) {return read_error(str);}
		}
		else {
			cout << "Unrecognized city keyword in input file: " << str << endl;
			return 0;
		}
		return 1;
	}
	vector3d get_car_size() const {return CAR_SIZE*road_width;}
}; // city_params_t

city_params_t city_params;


vector3d const get_query_xlate() {
	return vector3d((world_mode == WMODE_INF_TERRAIN) ? vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0) : zero_vector);
}

bool check_bcube_sphere_coll(cube_t const &bcube, point const &sc, float radius, bool xy_only) {
	return (xy_only ? sphere_cube_intersect_xy(sc, radius, bcube) : sphere_cube_intersect(sc, radius, bcube));
}
template<typename T> bool check_bcubes_sphere_coll(vector<T> const &bcubes, point const &sc, float radius, bool xy_only) {
	for (auto i = bcubes.begin(); i != bcubes.end(); ++i) {
		if (check_bcube_sphere_coll(*i, sc, radius, xy_only)) return 1;
	}
	return 0;
}

template<typename T> void get_all_bcubes(vector<T> const &v, vector<cube_t> &bcubes) {
	for (auto i = v.begin(); i != v.end(); ++i) {bcubes.push_back(*i);}
}

template<typename T> static void add_flat_road_quad(T const &r, quad_batch_draw &qbd, colorRGBA const &color, float ar) { // z1 == z2
	float const z(r.z1());
	point const pts[4] = {point(r.x1(), r.y1(), z), point(r.x2(), r.y1(), z), point(r.x2(), r.y2(), z), point(r.x1(), r.y2(), z)};
	qbd.add_quad_pts(pts, color, plus_z, r.get_tex_range(ar));
}

float smooth_interp(float a, float b, float mix) {
	mix = mix * mix * (3.0 - 2.0 * mix); // cubic Hermite interoplation (smoothstep)
	return mix*a + (1.0 - mix)*b;
}

bool is_isect(unsigned type) {return (type >= TYPE_ISEC2 && type <= TYPE_ISEC4);}
int encode_neg_ix(unsigned ix) {return -(int(ix)+1);}
unsigned decode_neg_ix(int ix) {assert(ix < 0); return -(ix+1);}

class road_mat_mgr_t {

	bool inited;
	unsigned tids[NUM_RD_TIDS], sl_tid;

public:
	road_mat_mgr_t() : inited(0), sl_tid(0) {}

	void ensure_road_textures() {
		if (inited) return;
		timer_t timer("Load Road Textures");
		string const img_names[NUM_RD_TIDS] = {"sidewalk.jpg", "straight_road.jpg", "bend_90.jpg", "int_3_way.jpg", "int_4_way.jpg", "parking_lot.png"};
		float const aniso[NUM_RD_TIDS] = {4.0, 16.0, 8.0, 8.0, 8.0, 4.0};
		for (unsigned i = 0; i < NUM_RD_TIDS; ++i) {tids[i] = get_texture_by_name(("roads/" + img_names[i]), 0, 0, 1, aniso[i]);}
		sl_tid = get_texture_by_name("roads/traffic_light.png");
		inited = 1;
	}
	void set_texture(unsigned type) {
		assert(type < NUM_RD_TYPES);
		ensure_road_textures();
		select_texture(tids[type]);
	}
	void set_stoplight_texture() {
		ensure_road_textures();
		select_texture(sl_tid);
	}
};

road_mat_mgr_t road_mat_mgr;

void city_shader_setup(shader_t &s, cube_t const &lights_bcube, bool use_dlights, bool use_smap, int use_bmap) {

	use_dlights &= !lights_bcube.is_zero_area();
	setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, use_dlights, 0, 0, use_smap, use_bmap, 0, use_dlights, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1

	if (use_dlights) {
		s.add_uniform_vector3d("scene_llc",   lights_bcube.get_llc()); // reset with correct values
		s.add_uniform_vector3d("scene_scale", lights_bcube.get_size());
		s.add_uniform_float("LT_DIR_FALLOFF", CITY_LIGHT_FALLOFF); // smooth falloff for car headlights and streetlights
	}
	if (use_smap) {
		s.add_uniform_float("z_bias", cobj_z_bias);
		s.add_uniform_float("pcf_offset", 8.0*shadow_map_pcf_offset);
		s.add_uniform_float("dlight_pcf_offset", 0.0005);
	}
}

struct draw_state_t {
	shader_t s;
	vector3d xlate;
protected:
	bool use_smap, use_bmap, shadow_only, use_dlights, emit_now;
	cube_t lights_bcube;
	point_sprite_drawer_sized light_psd; // for car/traffic lights
public:
	draw_state_t() : xlate(zero_vector), use_smap(0), use_bmap(0), shadow_only(0), use_dlights(0), emit_now(0) {}
	virtual void draw_unshadowed() {}
	void begin_tile(point const &pos) {emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));}

	void pre_draw(vector3d const &xlate_, cube_t const &lights_bcube_, bool use_dlights_, bool shadow_only_, bool always_setup_shader) {
		xlate        = xlate_;
		shadow_only  = shadow_only_;
		lights_bcube = lights_bcube_;
		use_dlights  = (use_dlights_ && !shadow_only);
		use_smap     = (shadow_map_enabled() && !shadow_only);
		if (!use_smap && !always_setup_shader) return;
		city_shader_setup(s, lights_bcube, use_dlights, use_smap, (use_bmap && !shadow_only));
	}
	virtual void post_draw() {
		emit_now = 0;
		if (use_smap) {s.end_shader();}
		if (shadow_only) {s.begin_color_only_shader();} else {city_shader_setup(s, lights_bcube, use_dlights, 0, use_bmap);} // no smap
		draw_unshadowed();
		s.end_shader();
	}
	void draw_and_clear_light_flares() {
		if (light_psd.empty()) return; // no lights to draw
		enable_blend();
		set_additive_blend_mode();
		light_psd.draw_and_clear(BLUR_TEX, 0.0, 0, 1, 0.005); // use geometry shader for unlimited point size
		set_std_blend_mode();
		disable_blend();
	}
	bool check_cube_visible(cube_t const &bc, float dist_scale=1.0, bool shadow_only=0) const {
		cube_t const bcx(bc + xlate);

		if (dist_scale > 0.0) {
			float const dmax(shadow_only ? camera_pdu.far_ : dist_scale*get_draw_tile_dist());
			if (!dist_less_than(camera_pdu.pos, bcx.closest_pt(camera_pdu.pos), dmax)) return 0;
		}
		return camera_pdu.cube_visible(bcx);
	}
	static void set_cube_pts(cube_t const &c, float z1, float z2, bool d, bool D, point p[8]) {
		p[0][!d] = p[4][!d] = c.d[!d][1]; p[0][d] = p[4][d] = c.d[d][ D]; p[0].z = z1; p[4].z = z2; // front right
		p[1][!d] = p[5][!d] = c.d[!d][0]; p[1][d] = p[5][d] = c.d[d][ D]; p[1].z = z1; p[5].z = z2; // front left
		p[2][!d] = p[6][!d] = c.d[!d][0]; p[2][d] = p[6][d] = c.d[d][!D]; p[2].z = z1; p[6].z = z2; // back  left
		p[3][!d] = p[7][!d] = c.d[!d][1]; p[3][d] = p[7][d] = c.d[d][!D]; p[3].z = z1; p[7].z = z2; // back  right
	}
	static void rotate_pts(point const &center, float sine_val, float cos_val, int d, int e, point p[8]) {
		for (unsigned i = 0; i < 8; ++i) {
			point &v(p[i]); // rotate p[i]
			v -= center; // translate to origin
			float const a(v[d]*cos_val - v[e]*sine_val), b(v[e]*cos_val + v[d]*sine_val);
			v[d] = a; v[e] = b; // rotate
			v += center; // translate back
		}
	}
	void draw_cube(quad_batch_draw &qbd, bool d, bool D, color_wrapper const &cw, point const &center, point const p[8]) const {
		vector3d const cview_dir((camera_pdu.pos - xlate) - center);
		float const sign((d^D) ? -1.0 : 1.0);
		vector3d const top_n  (cross_product((p[2] - p[1]), (p[0] - p[1]))*sign); // Note: normalization not needed
		vector3d const front_n(cross_product((p[5] - p[1]), (p[0] - p[1]))*sign);
		vector3d const right_n(cross_product((p[6] - p[2]), (p[1] - p[2]))*sign);
		if (dot_product(cview_dir, top_n) > 0) {qbd.add_quad_pts(p+4, cw,  top_n);} // top
		//else                                   {qbd.add_quad_pts(p+0, cw, -top_n);} // bottom - not actually drawn
		if (dot_product(cview_dir, front_n) > 0) {point const pts[4] = {p[0], p[1], p[5], p[4]}; qbd.add_quad_pts(pts, cw,  front_n);} // front
		else                                     {point const pts[4] = {p[2], p[3], p[7], p[6]}; qbd.add_quad_pts(pts, cw, -front_n);} // back
		if (dot_product(cview_dir, right_n) > 0) {point const pts[4] = {p[1], p[2], p[6], p[5]}; qbd.add_quad_pts(pts, cw,  right_n);} // right
		else                                     {point const pts[4] = {p[3], p[0], p[4], p[7]}; qbd.add_quad_pts(pts, cw, -right_n);} // left
	}
	bool add_light_flare(point const &flare_pos, vector3d const &n, colorRGBA const &color, float alpha, float radius) {
		point pos(xlate + flare_pos);
		vector3d const view_dir((camera_pdu.pos - pos).get_norm());
		float dp(dot_product(n, view_dir));
		pos += 0.75*radius*view_dir; // move toward the camera, away from the stoplight, to prevent clipping
		if (dp < 0.05) return 0; // back facing, skip
		light_psd.add_pt(sized_vert_t<vert_color>(vert_color(pos, colorRGBA(color, dp*alpha)), radius));
		return 1;
	}
}; // draw_state_t

class city_road_gen_t;

struct car_t {
	cube_t bcube, prev_bcube;
	bool dim, dir, stopped_at_light, entering_city;
	unsigned char cur_road_type, color_id, turn_dir, front_car_turn_dir, model_id;
	unsigned short cur_city, cur_road, cur_seg;
	float height, dz, rot_z, turn_val, cur_speed, max_speed;

	car_t() : bcube(all_zeros), dim(0), dir(0), stopped_at_light(0), entering_city(0), cur_road_type(0), color_id(0), turn_dir(TURN_NONE), front_car_turn_dir(TURN_UNSPEC),
		model_id(0), cur_city(0), cur_road(0), cur_seg(0), height(0.0), dz(0.0), rot_z(0.0), turn_val(0.0), cur_speed(0.0), max_speed(0.0) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	point get_center() const {return bcube.get_cube_center();}
	unsigned get_orient() const {return (2*dim + dir);}
	float get_max_speed() const {return ((cur_city == CONN_CITY_IX) ? CONN_ROAD_SPEED_MULT : 1.0)*max_speed;}
	float get_length() const {return (bcube.d[ dim][1] - bcube.d[ dim][0]);}
	float get_width () const {return (bcube.d[!dim][1] - bcube.d[!dim][0]);}
	bool is_almost_stopped() const {return (cur_speed < 0.1*max_speed);}
	bool is_stopped () const {return (cur_speed == 0.0);}
	bool is_parked  () const {return (max_speed == 0.0);}
	bool in_isect   () const {return is_isect(cur_road_type);}
	void park() {cur_speed = max_speed = 0.0;}
	float get_turn_rot_z(float dist_to_turn) const {return (1.0 - CLIP_TO_01(4.0f*fabs(dist_to_turn)/city_params.road_width));}

	string str() const {
		std::ostringstream oss;
		oss << "Car " << TXT(dim) << TXT(dir) << TXT(cur_city) << TXT(cur_road) << TXT(cur_seg) << TXT(dz) << TXT(max_speed) << TXT(cur_speed)
			<< "cur_road_type=" << unsigned(cur_road_type) << " color=" << unsigned(color_id) << " bcube=" << bcube.str();
		return oss.str();
	}
	void move(float speed_mult) {
		prev_bcube = bcube;
		if (is_stopped()) return;
		assert(speed_mult >= 0.0 && cur_speed > 0.0 && cur_speed <= CONN_ROAD_SPEED_MULT*max_speed); // Note: must be valid for connector road => city transitions
		float dist(cur_speed*speed_mult);
		if (dz != 0.0) {dist *= min(1.25, max(0.75, (1.0 - 0.5*dz/get_length())));} // slightly faster down hills, slightly slower up hills
		min_eq(dist, 0.25f*city_params.road_width); // limit to half a car length to prevent cars from crossing an intersection in a single frame
		move_by((dir ? 1.0 : -1.0)*dist);
	}
	void accelerate(float mult=0.02) {cur_speed = min(get_max_speed(), (cur_speed + mult*fticks*max_speed));}
	void decelerate(float mult=0.05) {cur_speed = max(0.0f, (cur_speed - mult*fticks*max_speed));}
	void decelerate_fast() {decelerate(10.0);} // Note: large decel to avoid stopping in an intersection
	void stop() {cur_speed = 0.0;} // immediate stop
	void move_by(float val) {bcube.d[dim][0] += val; bcube.d[dim][1] += val;}
	bool check_collision(car_t &c, city_road_gen_t const &road_gen);
};

struct comp_car_road_then_pos {
	vector3d const &xlate;
	comp_car_road_then_pos(vector3d const &xlate_) : xlate(xlate_) {}

	bool operator()(car_t const &c1, car_t const &c2) const { // sort spatially for collision detection and drawing
		if (c1.cur_city != c2.cur_city) return (c1.cur_city < c2.cur_city);
		if (c1.cur_road != c2.cur_road) return (c1.cur_road < c2.cur_road);
		if (c1.is_parked() != c2.is_parked()) {return c2.is_parked();} // parked cars last
		
		if (c1.is_parked()) { // sort parked cars back to front relative to camera so that alpha blending works
			return (p2p_dist_sq((c1.bcube.get_cube_center() + xlate), camera_pdu.pos) > p2p_dist_sq((c2.bcube.get_cube_center() + xlate), camera_pdu.pos));
		}
		return (c1.bcube.d[c1.dim][c1.dir] < c2.bcube.d[c2.dim][c2.dir]); // compare front end of car (used for collisions)
	}
};

struct rect_t {
	unsigned x1, y1, x2, y2;
	rect_t() : x1(0), y1(0), x2(0), y2(0) {}
	rect_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
	bool is_valid() const {return (x1 < x2 && y1 < y2);}
	unsigned get_area() const {return (x2 - x1)*(y2 - y1);}
	bool operator== (rect_t const &r) const {return (x1 == r.x1 && y1 == r.y1 && x2 == r.x2 && y2 == r.y2);}
	bool has_overlap(rect_t const &r) const {return (x1 < r.x2 && y1 < r.y2 && r.x1 < x2 && r.y1 < y2);}
};
struct flatten_op_t : public rect_t {
	float z1, z2;
	bool dim;
	unsigned border;
	flatten_op_t() : z1(0.0), z2(0.0), dim(0), border(0) {}
	flatten_op_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_, float z1_, float z2_, bool dim_, unsigned border_) :
		rect_t(x1_, y1_, x2_, y2_), z1(z1_), z2(z2_), dim(dim_), border(border_) {}
};

struct road_t : public cube_t {
	unsigned road_ix;
	bool dim; // dim the road runs in
	bool slope; // 0: z1 applies to first (lower) point; 1: z1 applies to second (upper) point

	road_t(cube_t const &c, bool dim_, bool slope_=0, unsigned road_ix_=0) : cube_t(c), road_ix(road_ix_), dim(dim_), slope(slope_) {}
	road_t(point const &s, point const &e, float width, bool dim_, bool slope_=0, unsigned road_ix_=0) : road_ix(road_ix_), dim(dim_), slope(slope_) {
		assert(s != e);
		assert(width > 0.0);
		vector3d const dw(0.5*width*cross_product((e - s), plus_z).get_norm());
		point const pts[4] = {(s - dw), (s + dw), (e + dw), (e - dw)};
		set_from_points(pts, 4);
	}
	float get_length   () const {return (d[dim][1] - d[dim][0]);}
	float get_slope_val() const {return get_dz()/get_length();}
};

struct road_seg_t : public road_t {
	unsigned short road_ix, conn_ix[2], conn_type[2]; // {dim=0, dim=1}
	void init_ixs() {conn_ix[0] = conn_ix[1] = 0; conn_type[0] = conn_type[1] = CONN_TYPE_NONE;}
	road_seg_t(road_t const &r, unsigned rix) : road_t(r), road_ix(rix) {init_ixs();}
	road_seg_t(cube_t const &c, unsigned rix, bool dim_, bool slope_=0) : road_t(c, dim_, slope_), road_ix(rix) {init_ixs();}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, -ar, (dim ? -1.0 : 1.0), 0, dim);}

	void add_road_quad(quad_batch_draw &qbd, colorRGBA const &color, float ar) const { // specialized here for sloped roads
		if (z1() == z2()) {add_flat_road_quad(*this, qbd, color, ar); return;}
		bool const s(slope ^ dim);
		point pts[4] = {point(x1(), y1(), d[2][!s]), point(x2(), y1(), d[2][!s]), point(x2(), y2(), d[2][ s]), point(x1(), y2(), d[2][ s])};
		if (!dim) {swap(pts[0].z, pts[2].z);}
		vector3d const normal(cross_product((pts[2] - pts[1]), (pts[0] - pts[1])).get_norm());
		qbd.add_quad_pts(pts, color, normal, get_tex_range(ar));
	}
};

namespace stoplight_ns {

	enum {GREEN_LIGHT=0, YELLOW_LIGHT, RED_LIGHT}; // colors, unused (only have stop and go states anyway)
	enum {EGL=0, EGWG, WGL, NGL, NGSG, SGL, NUM_STATE}; // E=car moving east, W=west, N=sorth, S=south, G=straight|right, L=left turn
	float const state_times[NUM_STATE] = {5.0, 6.0, 5.0, 5.0, 6.0, 5.0}; // in seconds
	unsigned const st_r_orient_masks[NUM_STATE] = {2, 3, 1, 8, 12, 4}; // {W=1, E=2, S=4, N=8}, for straight and right turns
	unsigned const left_orient_masks[NUM_STATE] = {2, 0, 1, 8, 0,  4}; // {W=1, E=2, S=4, N=8}, for left turns only
	unsigned const to_right  [4] = {3, 2, 0, 1}; // {N, S, W, E}
	unsigned const to_left   [4] = {2, 3, 1, 0}; // {S, N, E, W}
	unsigned const other_lane[4] = {1, 0, 3, 2}; // {E, W, N, S}

	rand_gen_t stoplight_rgen;

	class stoplight_t {
		uint8_t num_conn, conn, cur_state;
		bool at_conn_road; // longer light times in this case
		mutable bool blocked[4]; // Note: 4 bit flags corresponding to conn bits; mutable because it's set during car update logic, where roads are supposed to be const
		float cur_state_ticks;

		void next_state() {
			++cur_state;
			if (cur_state == NUM_STATE) {cur_state = 0;} // wraparound
		}
		void advance_state() {
			if (num_conn == 4) {next_state();} // all states are valid for 4-way intersections
			else { // 3-way intersection
				assert(num_conn == 3);
				while (1) {
					next_state();
					bool valid(0);
					switch (conn) {
					case 7 : {bool const allow[6] = {0,1,1,1,0,0}; valid = allow[cur_state]; break;} // no +y / S
					case 11: {bool const allow[6] = {1,1,0,0,0,1}; valid = allow[cur_state]; break;} // no -y / N
					case 13: {bool const allow[6] = {1,0,0,1,1,0}; valid = allow[cur_state]; break;} // no +x / W
					case 14: {bool const allow[6] = {0,0,1,0,1,1}; valid = allow[cur_state]; break;} // no -x / E
					default: assert(0);
					}
					if (valid) break;
				} // end while
			}
			cur_state_ticks = 0.0; // reset for this state
		}
		void run_update_logic() {
			assert(cur_state < NUM_STATE);
			if (cur_state_ticks > get_cur_state_time_secs()) {advance_state();} // time to update to next state
		}
		float get_cur_state_time_secs() const {return (at_conn_road ? 2.0 : 1.0)*TICKS_PER_SECOND*state_times[cur_state];}
	public:
		stoplight_t(bool at_conn_road_) : num_conn(0), conn(0), cur_state(RED_LIGHT), at_conn_road(at_conn_road_), cur_state_ticks(0.0) {UNROLL_4X(blocked[i_] = 0;)}
		void mark_blocked(bool dim, bool dir) const {blocked[2*dim + dir] = 1;} // Note: not actually const, but blocked is mutable
		bool is_blocked(bool dim, bool dir) const {return (blocked[2*dim + dir] != 0);}

		void init(uint8_t num_conn_, uint8_t conn_) {
			num_conn = num_conn_; conn = conn_;
			if (num_conn == 2) return; // nothing else to do
			cur_state = stoplight_rgen.rand() % NUM_STATE; // start at a random state
			advance_state(); // make sure cur_state is valid
			cur_state_ticks = get_cur_state_time_secs()*stoplight_rgen.rand_float(); // start at a random time within this state
		}
		void next_frame() {
			UNROLL_4X(blocked[i_] = 0;)
			if (num_conn == 2) return; // nothing else to do
			cur_state_ticks += fticks;
			run_update_logic();
		}
		bool red_light(bool dim, bool dir, unsigned turn) const {
			assert(cur_state < NUM_STATE);
			assert(turn == TURN_LEFT || turn == TURN_RIGHT || turn == TURN_NONE);
			if (num_conn == 2) return 0; // 2-way intersection, no cross traffic
			unsigned const mask(((turn == TURN_LEFT) ? left_orient_masks : st_r_orient_masks)[cur_state]);
			unsigned const orient(2*dim + dir); // {W, E, S, N}
			return (((1<<orient) & mask) == 0);
		}
		unsigned get_light_state(bool dim, bool dir, unsigned turn) const { // 0=green, 1=yellow, 2=red
			if (red_light(dim, dir, turn)) return RED_LIGHT;
			if (num_conn == 2) return GREEN_LIGHT;
			stoplight_t future_self(*this);
			float const yellow_light_time(2.0);
			future_self.cur_state_ticks += TICKS_PER_SECOND*yellow_light_time;
			future_self.run_update_logic();
			return (future_self.red_light(dim, dir, turn) ? YELLOW_LIGHT : GREEN_LIGHT);
		}
		bool check_int_clear(car_t const &car) const { // check for cars on other lanes blocking the intersection
			unsigned const orient(car.get_orient());  // {W, E, S, N}
			switch (car.turn_dir) {
			case TURN_NONE:  return (!blocked[to_right[orient]] && !blocked[to_left[orient]]); // straight
			case TURN_LEFT:  return (!blocked[to_right[orient]] && !blocked[to_left[orient]] && !blocked[other_lane[orient]]);
			case TURN_RIGHT: return (!blocked[to_right[orient]]);
			}
			return 1;
		}
		bool can_turn_right_on_red(car_t const &car) const { // check for legal right on red (no other lanes turning into the road to our right)
			if (car.turn_dir != TURN_RIGHT) return 0;
			unsigned const orient(car.get_orient());  // {W, E, S, N}
			if (!red_light(!car.dim, (to_right  [orient] & 1), TURN_NONE)) return 0; // traffic to our left has a green or yellow light for going straight
			if (!red_light( car.dim, (other_lane[orient] & 1), TURN_LEFT)) return 0; // opposing traffic has a green or yellow light for turning left
			// Note: there are no U-turns, so we don't need to worry about 
			return 1; // can turn right on red to_left[orient]
		}
		colorRGBA get_stoplight_color(bool dim, bool dir, unsigned turn) const {return stoplight_colors[get_light_state(dim, dir, turn)];}
	};
} // end stoplight_ns

struct road_isec_t : public cube_t {
	uint8_t num_conn, conn; // connected roads in {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
	short rix_xy[2], conn_ix[4]; // pos=cur city road, neg=global road; always segment ix
	stoplight_ns::stoplight_t stoplight; // Note: not always needed, maybe should be by pointer/index?

	road_isec_t(cube_t const &c, int rx, int ry, uint8_t conn_, bool at_conn_road) : cube_t(c), conn(conn_), stoplight(at_conn_road) {
		rix_xy[0] = rx; rix_xy[1] = ry; conn_ix[0] = conn_ix[1] = conn_ix[2] = conn_ix[3] = 0;
		if (conn == 15) {num_conn = 4;} // 4-way
		else if (conn == 7 || conn == 11 || conn == 13 || conn == 14) {num_conn = 3;} // 3-way
		else if (conn == 5 || conn == 6  || conn == 9  || conn == 10) {num_conn = 2;} // 2-way
		else {assert(0);}
		stoplight.init(num_conn, conn);
	}
	tex_range_t get_tex_range(float ar) const {
		switch (conn) {
		case 5 : return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 0); // 2-way: MX
		case 6 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 2-way: R0
		case 9 : return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 2-way: MXMY
		case 10: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 0); // 2-way: MY
		case 7 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 3-way: R0
		case 11: return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 3-way: MY
		case 13: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 1); // 3-way: R90MY
		case 14: return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 1); // 3-way: R90MX
		case 15: return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 4-way: R0
		default: assert(0);
		}
		return tex_range_t(0.0, 0.0, 1.0, 1.0); // never gets here
	}
	void make_4way() {num_conn = 4; conn = 15;}
	void next_frame() {stoplight.next_frame();}
	bool is_global_conn_int() const {return (rix_xy[0] < 0 || rix_xy[1] < 0);}
	bool red_light(car_t const &car) const {return stoplight.red_light(car.dim, car.dir, car.turn_dir);}
	bool red_or_yellow_light(car_t const &car) const {return (stoplight.get_light_state(car.dim, car.dir, car.turn_dir) != stoplight_ns::GREEN_LIGHT);}

	bool can_go_now(car_t const &car) const {
		if (!stoplight.check_int_clear(car)) return 0; // intersection not clear
		if (!red_or_yellow_light(car)) return 1; // green light
		return stoplight.can_turn_right_on_red(car); // stopped at light
	}
	bool has_left_turn_signal(unsigned orient) const {
		if (num_conn == 2) return 0; // never
		if (num_conn == 4) return 1; // always
		assert(num_conn == 3);
		switch (conn) {
		case 7 : return (orient == 1 || orient == 2);
		case 11: return (orient == 0 || orient == 3);
		case 13: return (orient == 0 || orient == 2);
		case 14: return (orient == 1 || orient == 3);
		default: assert(0);
		}
		return 0;
	}

	void draw_sl_block(quad_batch_draw &qbd, draw_state_t &dstate, point p[4], float h, unsigned state, bool draw_unlit, float flare_alpha, vector3d const &n, tex_range_t const &tr) const {
		for (unsigned j = 0; j < 3; ++j) {
			colorRGBA const &color(stoplight_colors[j]);

			if (j == state) {
				qbd.add_quad_pts(p, color, n, tr);
				if (flare_alpha > 0.0) {dstate.add_light_flare(0.25*(p[0] + p[1] + p[2] + p[3]), n, color, flare_alpha, 2.0*h);}
			}
			else if (draw_unlit) {
				qbd.add_quad_pts(p, (color + WHITE)*0.05, n, tr);
			}
			for (unsigned e = 0; e < 4; ++e) {p[e].z += 1.2*h;}
		}
	}
	void draw_stoplights(quad_batch_draw &qbd, draw_state_t &dstate, bool shadow_only) const {
		if (num_conn == 2) return; // no stoplights
		if (!dstate.check_cube_visible(*this, 0.16, shadow_only)) return; // dist_scale=0.12
		point const center(get_cube_center() + dstate.xlate);
		float const dist_val(shadow_only ? 0.0 : p2p_dist(camera_pdu.pos, center)/get_draw_tile_dist());
		vector3d const cview_dir(camera_pdu.pos - center);
		float const sz(0.03*city_params.road_width), h(1.0*sz);
		color_wrapper cw; cw.set_c4(BLACK);

		for (unsigned n = 0; n < 4; ++n) { // {-x, +x, -y, +y} = {W, E, S, N} facing = car traveling {E, W, N, S}
			if (!(conn & (1<<n))) continue; // no road in this dir
			bool const dim((n>>1) != 0), dir((n&1) == 0), side((dir^dim^1) != 0); // Note: dir is inverted here to represent car dir
			float const zbot(z1() + 2.0*h);
			float const pos(d[dim][!dir] + (dir ? sz : -sz)); // location in road dim
			float const v1(d[!dim][side]), v2(v1 + (side ? -sz : sz)); // location in other dim
			// draw base
			unsigned const num_segs(has_left_turn_signal(n) ? 6 : 3);
			float const sl_top(zbot + 1.2*h*num_segs), sl_lo(min(v1, v2) - 0.25*sz), sl_hi(max(v1, v2) + 0.25*sz);

			if (dist_val > 0.06) { // draw front face only
				point pts[4];
				pts[0][dim]  = pts[1][dim]  = pts[2][dim] = pts[3][dim] = pos;
				pts[0][!dim] = pts[3][!dim] = sl_lo;
				pts[1][!dim] = pts[2][!dim] = sl_hi;
				pts[0].z = pts[1].z = z1();
				pts[2].z = pts[3].z = sl_top;
				qbd.add_quad_pts(pts, cw,  (dim ? (dir ? plus_y : -plus_y) : (dir ? plus_x : -plus_x))); // Note: normal doesn't really matter since color is black
			}
			else {
				cube_t c;
				c.z1() = z1(); c.z2() = sl_top;
				c.d[ dim][0] = pos - (dir ? -0.04 : 0.5)*sz; c.d[dim][1] = pos + (dir ? 0.5 : -0.04)*sz;
				c.d[!dim][0] = sl_lo; c.d[!dim][1] = sl_hi;
				point pts[8];
				dstate.set_cube_pts(c, c.z1(), c.z2(), dim, dir, pts);
				dstate.draw_cube(qbd, dim, dir, cw, c.get_cube_center(), pts); // Note: uses traffic light texture, but color is back so it's all black anyway
			}
			if (shadow_only)    continue; // no lights in shadow pass
			if (dist_val > 0.1) continue; // too far away
			vector3d normal(zero_vector);
			normal[dim] = (dir ? -1.0 : 1.0);
			if (dot_product(normal, cview_dir) < 0.0) continue; // back facing, don't draw the lights
			// draw straight/line turn light
			point p[4];
			p[0][dim] = p[1][dim] = p[2][dim] = p[3][dim] = pos;
			p[0][!dim] = p[3][!dim] = v1; p[1][!dim] = p[2][!dim] = v2;
			p[0].z = p[1].z = zbot; p[2].z = p[3].z = zbot + h;
			bool const draw_detail(dist_val < 0.05); // only when very close
			draw_sl_block(qbd, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_NONE), draw_detail, draw_detail, normal, tex_range_t(0.0, 0.0, 0.5, 1.0));
			// draw left turn light (upper light section)
			if (has_left_turn_signal(n)) {
				draw_sl_block(qbd, dstate, p, h, stoplight.get_light_state(dim, dir, TURN_LEFT), draw_detail, 0.5*draw_detail, normal, tex_range_t(1.0, 0.0, 0.5, 1.0));
			}
		} // for n
	}
};

struct road_plot_t : public cube_t {
	road_plot_t(cube_t const &c) : cube_t(c) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, ar, ar);}
};

struct parking_lot_t : public cube_t {
	bool dim, dir;
	unsigned short row_sz, num_rows;
	parking_lot_t(cube_t const &c, bool dim_, bool dir_, unsigned row_sz_=0, unsigned num_rows_=0) : cube_t(c), dim(dim_), dir(dir_), row_sz(row_sz_), num_rows(num_rows_) {}
	
	tex_range_t get_tex_range(float ar) const {
		bool const d(!dim); // Note: R90
		float const xscale(1.0/(2.0*PARK_SPACE_WIDTH *city_params.get_car_size().y));
		float const yscale(1.0/(1.0*PARK_SPACE_LENGTH*city_params.get_car_size().x));
		float const dx(get_dx()), dy(get_dy()), tx(0.24), ty(0.0); // x=cols, y=rows
		return tex_range_t(tx, ty, (xscale*(d ? dy : dx) + tx), (yscale*(d ? dx : dy) + ty), 0, d);
	}
};

bool is_night() {return (light_factor < 0.5);} // for car headlights and streetlights

namespace streetlight_ns {

	colorRGBA const pole_color(BLACK); // so we don't have to worry about shadows
	colorRGBA const light_color(1.0, 0.9, 0.7, 1.0);
	float const light_height = 0.5; // in units of road width
	float const pole_radius  = 0.015;
	float const light_radius = 0.025;
	float const light_dist   = 3.0;

	struct streetlight_t {
		point pos; // bottom point
		vector3d dir;

		streetlight_t(point const &pos_, vector3d const &dir_) : pos(pos_), dir(dir_) {}

		point get_lpos() const {
			float const height(light_height*city_params.road_width);
			return (pos + vector3d(0.0, 0.0, 1.1*height) + 0.4*height*dir);
		}
		void draw(shader_t &s, vector3d const &xlate, bool shadow_only) const { // Note: translate has already been applied as a transform
			float const height(light_height*city_params.road_width);
			point const center(pos + xlate + vector3d(0.0, 0.0, 0.5*height));
			if (shadow_only && !dist_less_than(camera_pdu.pos, center, 0.8*camera_pdu.far_)) return;
			if (!camera_pdu.sphere_visible_test(center, height)) return; // VFC
			float const dist_val(shadow_only ? 0.06 : p2p_dist(camera_pdu.pos, center)/get_draw_tile_dist());
			if (dist_val > 0.2) return; // too far
			if (!s.is_setup()) {s.begin_color_only_shader();} // likely only needed for shadow pass, could do better with this
			float const pradius(pole_radius*city_params.road_width), lradius(light_radius*city_params.road_width);
			int const ndiv(max(4, min(N_SPHERE_DIV, int(0.5/dist_val))));
			point const top(pos + vector3d(0.0, 0.0, 0.96*height)), lpos(get_lpos()), arm_end(lpos + vector3d(0.0, 0.0, 0.025*height) - 0.06*height*dir);
			if (!shadow_only) {s.set_cur_color(pole_color);}
			draw_cylinder_at(pos, height, pradius, 0.7*pradius, min(ndiv, 24), 0); // vertical post, no ends
			if (dist_val < 0.12) {draw_fast_cylinder(top, arm_end, 0.5*pradius, 0.4*pradius, min(ndiv, 16), 0, 0);} // untextured, no ends
			
			if (shadow_only) {
				if (dist_less_than(camera_pdu.pos, (get_lpos() + xlate), 0.01*lradius)) return; // this is the light source, don't make it shadow itself
			}
			else {
				if (!is_night() && dist_val > 0.15) return; // too far
				if (is_night()) {s.set_color_e(light_color);} else {s.set_cur_color(light_color);} // emissive when lit
			}
			fgPushMatrix();
			translate_to(lpos);
			scale_by(lradius*vector3d(1.0+fabs(dir.x), 1.0+fabs(dir.y), 1.0)); // scale 2x in dir
			draw_sphere_vbo(all_zeros, 1.0, ndiv, 0); // untextured
			fgPopMatrix();
			if (!shadow_only && is_night()) {s.clear_color_e();}
		}
		void add_dlight(vector3d const &xlate, cube_t &lights_bcube) const {
			float const ldist(light_dist*city_params.road_width);
			if (!lights_bcube.contains_pt_xy(pos)) return; // not contained within the light volume
			float const height(light_height*city_params.road_width);
			point const lpos(get_lpos());
			if (!camera_pdu.sphere_visible_test((lpos + xlate), ldist)) return; // VFC
			min_eq(lights_bcube.z1(), (lpos.z - ldist));
			max_eq(lights_bcube.z2(), (lpos.z + ldist));
			dl_sources.push_back(light_source(ldist, lpos, lpos, light_color, 0, -plus_z, STREETLIGHT_BEAMWIDTH)); // points down
		}
	};
} // streetlight_ns

class heightmap_query_t {
protected:
	float *heightmap;
	unsigned xsize, ysize;

public:
	flatten_op_t last_flatten_op;

	heightmap_query_t() : heightmap(nullptr), xsize(0), ysize(0) {}
	heightmap_query_t(float *hmap, unsigned xsize_, unsigned ysize_) : heightmap(hmap), xsize(xsize_), ysize(ysize_) {}
	float get_x_value(int x) const {return get_xval(x - int(xsize)/2);} // convert from center to LLC
	float get_y_value(int y) const {return get_yval(y - int(ysize)/2);}
	int get_x_pos(float x) const {return (get_xpos(x) + int(xsize)/2);}
	int get_y_pos(float y) const {return (get_ypos(y) + int(ysize)/2);}
	float  get_height(unsigned x, unsigned y) const {return heightmap[y*xsize + x];} // Note: not bounds checked
	float &get_height(unsigned x, unsigned y)       {return heightmap[y*xsize + x];} // Note: not bounds checked
	bool is_normalized_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <  x2 && y1 <  y2 && x2 <= xsize && y2 <= ysize);}
	bool is_valid_region     (unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <= x2 && y1 <= y2 && x2 <= xsize && y2 <= ysize);}
	bool is_inside_terrain(int x, int y) const {return (x >= 0 && y >= 0 && x < (int)xsize && y < (int)ysize);}
	
	float get_height_at(float xval, float yval) const {
		int const x(get_x_pos(xval)), y(get_y_pos(yval));
		return (is_inside_terrain(x, y) ? get_height(x, y) : OUTSIDE_TERRAIN_HEIGHT);
	}
	float any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2, bool check_border=0) const {
		min_eq(x2, xsize); min_eq(y2, ysize); // clamp upper bound
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (check_border && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				if (get_height(x, y) < water_plane_z) return 1;
			}
		}
		return 0;
	}
	void flatten_region_to(cube_t const c, unsigned slope_width, bool decrease_only=0) {
		flatten_region_to(get_x_pos(c.x1()), get_y_pos(c.y1()), get_x_pos(c.x2()), get_y_pos(c.y2()), slope_width, (c.z1() - ROAD_HEIGHT), decrease_only);
	}
	void flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation, bool decrease_only=0) {
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = max((int)y1-(int)slope_width, 0); y < min(y2+slope_width, ysize); ++y) {
			for (unsigned x = max((int)x1-(int)slope_width, 0); x < min(x2+slope_width, xsize); ++x) {
				float &h(get_height(x, y));
				if (decrease_only && h < elevation) continue; // don't increase

				if (slope_width > 0) {
					float const dx(max(0, max(((int)x1 - (int)x), ((int)x - (int)x2 + 1))));
					float const dy(max(0, max(((int)y1 - (int)y), ((int)y - (int)y2 + 1))));
					h = smooth_interp(h, elevation, min(1.0f, sqrt(dx*dx + dy*dy)/slope_width));
				} else {h = elevation;}
			} // for x
		} // for y
	}
	float flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, unsigned border, bool stats_only=0, bool decrease_only=0) {
		if (!stats_only) {last_flatten_op = flatten_op_t(x1, y1, x2, y2, z1, z2, dim, border);} // cache for later replay
		assert(is_valid_region(x1, y1, x2, y2));
		if (x1 == x2 || y1 == y2) return 0.0; // zero area
		float const run_len(dim ? (y2 - y1) : (x2 - x1)), denom(1.0f/max(run_len, 1.0f)), dz(z2 - z1), border_inv(1.0/border);
		int const pad(border + 1U); // pad an extra 1 texel to handle roads misaligned with the texture
		unsigned px1(x1), py1(y1), px2(x2), py2(y2);
		float tot_dz(0.0);
		
		if (dim) {
			px1 = max((int)x1-pad, 0);
			px2 = min(x2+pad, xsize);
			py1 = max((int)y1-1, 0); // pad by 1 in road dim as well to blend with edge of city
			py2 = min(y2+1, ysize);
		}
		else {
			py1 = max((int)y1-pad, 0);
			py2 = min(y2+pad, ysize);
			px1 = max((int)x1-1, 0);
			px2 = min(x2+1, xsize);
		}
		for (unsigned y = py1; y < py2; ++y) {
			for (unsigned x = px1; x < px2; ++x) {
				float const t(((dim ? int(y - y1) : int(x - x1)) + ((dz < 0.0) ? 1 : -1))*denom); // bias toward the lower zval
				float const road_z(z1 + dz*t - ROAD_HEIGHT);
				float &h(get_height(x, y));
				if (decrease_only && h < road_z) continue; // don't increase
				float new_h;

				if (border > 0) {
					unsigned const dist(dim ? max(0, max(((int)x1 - (int)x - 1), ((int)x - (int)x2))) : max(0, max(((int)y1 - (int)y - 1), ((int)y - (int)y2))));
					new_h = smooth_interp(h, road_z, dist*border_inv);
				} else {new_h = road_z;}
				tot_dz += fabs(h - new_h);
				if (!stats_only) {h = new_h;} // apply the height change
			} // for x
		} // for y
		return tot_dz;
	}
	float flatten_for_road(road_t const &road, unsigned border, bool stats_only=0, bool decrease_only=0) {
		float const z_adj(ROAD_HEIGHT + 0.5*road.get_slope_val()*(road.dim ? DY_VAL : DX_VAL)); // account for a half texel of error for sloped roads
		unsigned const rx1(get_x_pos(road.x1())), ry1(get_y_pos(road.y1())), rx2(get_x_pos(road.x2())), ry2(get_y_pos(road.y2()));
		return flatten_sloped_region(rx1, ry1, rx2, ry2, road.d[2][road.slope]-z_adj, road.d[2][!road.slope]-z_adj, road.dim, border, stats_only, decrease_only);
	}
};

class city_plot_gen_t : public heightmap_query_t {

protected:
	int last_rgi;
	rand_gen_t rgen;
	vector<rect_t> used;
	vector<cube_t> plots; // same size as used
	cube_t bcube;

	bool overlaps_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		rect_t const cur(x1, y1, x2, y2);
		for (vector<rect_t>::const_iterator i = used.begin(); i != used.end(); ++i) {if (i->has_overlap(cur)) return 1;} // simple linear iteration
		return 0;
	}
	cube_t add_plot(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) {
		cube_t c;
		c.x1() = get_x_value(x1);
		c.x2() = get_x_value(x2);
		c.y1() = get_y_value(y1);
		c.y2() = get_y_value(y2);
		c.z1() = c.z2() = elevation;
		if (plots.empty()) {bcube = c;} else {bcube.union_with_cube(c);}
		plots.push_back(c);
		used.emplace_back(x1, y1, x2, y2);
		return c;
	}
	float get_avg_height(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		assert(is_normalized_region(x1, y1, x2, y2));
		float sum(0.0), denom(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				sum   += get_height(x, y);
				denom += 1.0;
			}
		}
		return sum/denom;
	}
	float get_rms_height_diff(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		float const avg(get_avg_height(x1, y1, x2, y2));
		float diff(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				float const delta(get_height(x, y) - avg);
				diff += delta*delta; // square the difference
			}
		}
		return diff;
	}
public:
	city_plot_gen_t() : last_rgi(0), bcube(all_zeros) {}

	void init(float *heightmap_, unsigned xsize_, unsigned ysize_) {
		heightmap = heightmap_; xsize = xsize_; ysize = ysize_;
		assert(heightmap != nullptr);
		assert(xsize > 0 && ysize > 0); // any size is okay
		if (rand_gen_index != last_rgi) {rgen.set_state(rand_gen_index, 12345); last_rgi = rand_gen_index;} // only when rand_gen_index changes
	}
	bool find_best_city_location(unsigned wmin, unsigned hmin, unsigned wmax, unsigned hmax, unsigned border, unsigned slope_width, unsigned num_samples,
		unsigned &cx1, unsigned &cy1, unsigned &cx2, unsigned &cy2)
	{
		assert(num_samples > 0);
		assert((wmax + 2*border) < xsize && (hmax + 2*border) < ysize); // otherwise the city can't fit in the map
		unsigned const num_iters(100*num_samples); // upper bound
		unsigned xend(xsize - wmax - 2*border + 1), yend(ysize - hmax - 2*border + 1); // max rect LLC, inclusive
		unsigned num_cands(0);
		float best_diff(0.0);

		for (unsigned n = 0; n < num_iters; ++n) { // find min RMS height change across N samples
			unsigned const x1(border + (rgen.rand()%xend)), y1(border + (rgen.rand()%yend));
			unsigned const x2(x1 + ((wmin == wmax) ? wmin : rgen.rand_int(wmin, wmax)));
			unsigned const y2(y1 + ((hmin == hmax) ? hmin : rgen.rand_int(hmin, hmax)));
			if (overlaps_used (x1-slope_width, y1-slope_width, x2+slope_width, y2+slope_width)) continue; // skip if plot expanded by slope_width overlaps an existing city
			if (any_underwater(x1, y1, x2, y2, CHECK_HEIGHT_BORDER_ONLY))   continue; // skip
			float const diff(get_rms_height_diff(x1, y1, x2, y2));
			if (num_cands == 0 || diff < best_diff) {cx1 = x1; cy1 = y1; cx2 = x2; cy2 = y2; best_diff = diff;}
			if (++num_cands == num_samples) break; // done
		} // for n
		if (num_cands == 0) return 0;
		//cout << "City cands: " << num_cands << ", diff: " << best_diff << ", loc: " << (cx1+cx2)/2 << "," << (cy1+cy2)/2 << endl;
		return 1; // success
	}
	float flatten_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float const *const height=nullptr) {
		float const elevation(height ? *height : get_avg_height(x1, y1, x2, y2));
		flatten_region_to(x1, y1, x2, y2, slope_width, elevation);
		return elevation;
	}
	bool check_plot_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (plots.empty()) return 0;
		point const query_pos(pos - get_query_xlate());
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
		return check_bcubes_sphere_coll(plots, query_pos, radius, xy_only);
	}
}; // city_plot_gen_t


class city_road_gen_t {

	struct range_pair_t {
		unsigned s, e; // Note: e is one past the end
		range_pair_t(unsigned s_=0, unsigned e_=0) : s(s_), e(e_) {}
		void update(unsigned v) {
			if (s == 0 && e == 0) {s = v;} // first insert
			else {assert(s < e && v >= e);} // v must strictly increase
			e = v+1; // one past the end
		}
	};

	class road_draw_state_t : public draw_state_t {
		quad_batch_draw qbd_batched[NUM_RD_TYPES], qbd_sl;
		float ar;

	public:
		road_draw_state_t() : ar(1.0) {}

		void pre_draw(vector3d const &xlate_, cube_t const &lights_bcube_, bool use_dlights_, bool shadow_only) {
			draw_state_t::pre_draw(xlate_, lights_bcube_, use_dlights_, shadow_only, 0); // always_setup_shader=0
			ar = city_params.get_road_ar();
		}
		virtual void draw_unshadowed() {
			for (unsigned i = 0; i < NUM_RD_TYPES; ++i) { // only unshadowed blocks
				road_mat_mgr.set_texture(i);
				qbd_batched[i].draw_and_clear();
			}
		}
		virtual void post_draw() {
			draw_state_t::post_draw();
			if (qbd_sl.empty()) return; // no stoplights to draw
			glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
			shader_t s;
			s.begin_simple_textured_shader(); // Note: no lighting
			road_mat_mgr.set_stoplight_texture();
			qbd_sl.draw_and_clear();
			s.end_shader();
			glDepthFunc(GL_LESS);
		}
		template<typename T> void add_road_quad(T const &r, quad_batch_draw &qbd, colorRGBA const &color) {add_flat_road_quad(r, qbd, color, ar);} // generic flat road case (plot/park)
		template<> void add_road_quad(road_seg_t  const &r, quad_batch_draw &qbd, colorRGBA const &color) {r.add_road_quad(qbd, color, ar);} // road segment
		
		template<typename T> void draw_road_region(vector<T> const &v, range_pair_t const &rp, quad_batch_draw &cache, unsigned type_ix) {
			assert(rp.s <= rp.e);
			assert(rp.e <= v.size());
			assert(type_ix < NUM_RD_TYPES);
			colorRGBA const color(road_colors[type_ix]);
			
			if (cache.empty()) { // generate and cache quads
				for (unsigned i = rp.s; i < rp.e; ++i) {add_road_quad(v[i], cache, color);}
			}
			if (emit_now) { // draw shadow blocks directly
				road_mat_mgr.set_texture(type_ix);
				cache.draw();
			} else {qbd_batched[type_ix].add_quads(cache);} // add non-shadow blocks for drawing later
		}
		void draw_stoplights(vector<road_isec_t> const &isecs, bool shadow_only) {
			for (auto i = isecs.begin(); i != isecs.end(); ++i) {i->draw_stoplights(qbd_sl, *this, shadow_only);}
		}
	}; // road_draw_state_t

	class parking_lot_manager_t {
	public: // road network needs access to parks for drawing
		vector<parking_lot_t> parks; // no, not really parks, but parking lots (the name "plots" was already taken)
	private:
		static bool has_bcube_int_xy(cube_t const &bcube, vector<cube_t> const &bcubes, float pad_dist=0.0) {
			cube_t tc(bcube);
			tc.expand_by(pad_dist);

			for (auto c = bcubes.begin(); c != bcubes.end(); ++c) {
				if (c->intersects_xy(tc)) return 1; // intersection
			}
			return 0;
		}
	public:
		void clear() {parks.clear();}

		void gen_parking(vector<road_plot_t> const &plots, vector<car_t> &cars, unsigned city_id) {
			if (city_params.min_park_spaces == 0 || city_params.min_park_rows == 0) return; // disable parking lots
			timer_t timer("Gen Parking Lots");
			vector3d const nom_car_size(city_params.get_car_size()); // {length, width, height}
			float const space_width(PARK_SPACE_WIDTH *nom_car_size.y); // add 50% extra space between cars
			float const space_len  (PARK_SPACE_LENGTH*nom_car_size.x); // space for car + gap for cars to drive through
			float const pad_dist   (1.0*nom_car_size.x); // one car length
			vector<cube_t> bcubes; // reused across calls
			rand_gen_t rgen;
			unsigned num_spaces(0), filled_spaces(0);
			parks.clear();
			rgen.set_state(city_id, 123);
			// cars
			car_t car;
			car.park();
			car.cur_city = city_id;
			car.cur_road_type = TYPE_PLOT;

			// generate 0-4 parking lots per plot, starting at the corners
			for (auto i = plots.begin(); i != plots.end(); ++i) {
				cube_t plot(*i);
				plot.expand_by_xy(-pad_dist);
				bcubes.clear();
				get_building_bcubes(plot, bcubes);
				if (bcubes.empty()) continue; // shouldn't happen, unless buildings are disabled; skip to avoid perf problems with an entire plot of parking lot
				unsigned const first_corner(rgen.rand()&3); // 0-3
				bool const car_dim(rgen.rand() & 1); // 0=cars face in X; 1=cars face in Y
				bool const car_dir(rgen.rand() & 1);
				float const xsz(car_dim ? space_width : space_len), ysz(car_dim ? space_len : space_width);
				//cout << "max_row_sz: " << floor(plot.get_size()[!car_dim]/space_width) << ", max_num_rows: " << floor(plot.get_size()[car_dim]/space_len) << endl;
				
				for (unsigned c = 0; c < 4; ++c) { // 4 corners, in random order
					unsigned const cix((first_corner + c) & 3), xdir(cix & 1), ydir(cix >> 1), wdir(car_dim ? xdir : ydir), rdir(car_dim ? ydir : xdir);
					float const dx(xdir ? -xsz : xsz), dy(ydir ? -ysz : ysz), dw(car_dim ? dx : dy), dr(car_dim ? dy : dx); // delta-wdith and delta-row
					point const corner_pos(plot.d[0][xdir], plot.d[1][ydir], (plot.z1() + 0.1*ROAD_HEIGHT)); // shift up slightly to avoid z-fighting
					assert(dw != 0.0 && dr != 0.0);
					parking_lot_t cand(cube_t(corner_pos, corner_pos), car_dim, car_dir, city_params.min_park_spaces, city_params.min_park_rows); // start as min size at the corner
					cand.d[!car_dim][!wdir] += cand.row_sz*dw;
					cand.d[ car_dim][!rdir] += cand.num_rows*dr;
					if (!plot.contains_cube_xy(cand)) {continue;} // can't fit a min size parking lot in this plot, so skip it (shouldn't happen)
					if (has_bcube_int_xy(cand, bcubes, pad_dist)) continue; // intersects a building - skip (can't fit min size parking lot)
					cand.z2() += plot.get_dz(); // probably unnecessary
					parking_lot_t park(cand);
					
					// try to add more parking spaces in a row
					for (; plot.contains_cube_xy(cand); ++cand.row_sz, cand.d[!car_dim][!wdir] += dw) {
						if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
						park = cand; // success: increase parking lot to this size
					}
					cand = park;
					// try to add more rows of parking spaces
					for (; plot.contains_cube_xy(cand); ++cand.num_rows, cand.d[car_dim][!rdir] += dr) {
						if (has_bcube_int_xy(cand, bcubes, pad_dist)) break; // intersects a building - done
						park = cand; // success: increase parking lot to this size
					}
					assert(park.row_sz >= city_params.min_park_spaces && park.num_rows >= city_params.min_park_rows);
					assert(park.get_dx() > 0.0 && park.get_dy() > 0.0);
					parks.push_back(park);
					//parks.back().expand_by_xy(0.5*pad_dist); // re-add half the padding for drawing (breaks texture coord alignment)
					bcubes.push_back(park); // add to list of blocker bcubes so that no later parking lots overlap this one
					num_spaces += park.row_sz*park.num_rows;

					// fill the parking lot with cars
					vector3d car_sz(nom_car_size);
					car.dim    = car_dim;
					car.dir    = car_dir;
					car.height = car_sz.z;
					if (car.dim) {swap(car_sz.x, car_sz.y);}
					point pos(corner_pos.x, corner_pos.y, (i->z2() + 0.5*car_sz.z));
					pos[ car_dim] += 0.5*dr + (car_dim ? 0.15 : -0.15)*fabs(dr); // offset for centerline, biased toward the front of the parking space
					float const car_density(rgen.rand_uniform(city_params.min_park_density, city_params.max_park_density));

					for (unsigned row = 0; row < park.num_rows; ++row) {
						pos[!car_dim] = corner_pos[!car_dim] + 0.5*dw; // half offset for centerline
						bool prev_was_bad(0);

						for (unsigned col = 0; col < park.row_sz; ++col) {
							if (prev_was_bad) {prev_was_bad = 0;} // previous car did a bad parking job, leave this space empty
							else if (rgen.rand_float() < car_density) { // only half the spaces are filled on average
								point cpos(pos);
								cpos[ car_dim] += 0.05*dr*rgen.rand_uniform(-1.0, 1.0); // randomness of front amount
								cpos[!car_dim] += 0.12*dw*rgen.rand_uniform(-1.0, 1.0); // randomness of side  amount
								
								if (col+1 != park.row_sz && (rgen.rand()&15) == 0) {// occasional bad parking job
									cpos[!car_dim] += dw*rgen.rand_uniform(0.3, 0.35);
									prev_was_bad = 1;
								} 
								car.bcube.set_from_point(cpos);
								car.bcube.expand_by(0.5*car_sz);
								cars.push_back(car);
								if ((rgen.rand()&7) == 0) {cars.back().dir ^= 1;} // pack backwards 1/8 of the time
								++filled_spaces;
							}
							pos[!car_dim] += dw;
						} // for col
						pos[car_dim] += dr;
					} // for row
					//cout << "plot: " << (i-plots.begin()) << ", b: " << bcubes.size() << ", dim: " << car_dim << ", dir: " << car_dir << ", row: " << park.row_sz << ", rows: " << park.num_rows << endl;
				} // for c
			} // for i
			cout << "parking lots: " << parks.size() << ", spaces: " << num_spaces << ", filled: " << filled_spaces << endl;
		}
	}; // parking_lot_manager_t

	class road_network_t {
		vector<road_t> roads; // full overlapping roads, for collisions, etc.
		vector<road_seg_t> segs; // non-overlapping road segments, for drawing with textures
		vector<road_isec_t> isecs[3]; // for drawing with textures: {4-way, 3-way, 2-way}
		vector<road_plot_t> plots; // plots of land that can hold buildings
		vector<streetlight_ns::streetlight_t> streetlights;
		parking_lot_manager_t parking_lot_mgr;
		cube_t bcube;
		vector<road_t> segments; // reused temporary
		set<unsigned> connected_to; // vector?
		map<uint64_t, unsigned> tile_to_block_map;
		unsigned city_id, cluster_id;
		//string city_name; // future work

		// use only for the global road network
		struct city_id_pair_t {
			unsigned id[2]; // lo, hi
			city_id_pair_t(unsigned c1, unsigned c2) {id[0] = c1; id[1] = c2;}
		};
		vector<city_id_pair_t> road_to_city; // indexed by road ID
		vector<vector<unsigned>> city_to_seg; // maps city_id to set of road segments connecting to that city

		static uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}

		struct cmp_by_tile { // not the most efficient solution, but no memory overhead
			bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
		};
		struct tile_block_t { // collection of road parts for a given tile
			range_pair_t ranges[NUM_RD_TYPES]; // {plot, seg, isec2, isec3, isec4, park_lot}
			quad_batch_draw quads[NUM_RD_TYPES];
			cube_t bcube;
			tile_block_t(cube_t const &bcube_) : bcube(bcube_) {}
		};
		vector<tile_block_t> tile_blocks;

		template<typename T> void add_tile_blocks(vector<T> &v, map<uint64_t, unsigned> &tile_to_block_map, unsigned type_ix) {
			assert(type_ix < NUM_RD_TYPES);
			sort(v.begin(), v.end(), cmp_by_tile());

			for (unsigned i = 0; i < v.size(); ++i) {
				uint64_t const tile_id(get_tile_id_for_cube(v[i]));
				auto it(tile_to_block_map.find(tile_id));
				unsigned block_id(0);
			
				if (it == tile_to_block_map.end()) { // not found, add new block
					tile_to_block_map[tile_id] = block_id = tile_blocks.size();
					tile_blocks.push_back(tile_block_t(v[i]));
				}
				else {block_id = it->second;}
				assert(block_id < tile_blocks.size());
				tile_blocks[block_id].ranges[type_ix].update(i);
				tile_blocks[block_id].bcube.union_with_cube(v[i]);
			} // for i
		}
	public:
		road_network_t() : bcube(all_zeros), city_id(CONN_CITY_IX), cluster_id(0) {} // global road network ctor
		road_network_t(cube_t const &bcube_, unsigned city_id_) : bcube(bcube_), city_id(city_id_), cluster_id(0) {bcube.d[2][1] += ROAD_HEIGHT;} // make it nonzero size
		cube_t const &get_bcube() const {return bcube;}
		void set_bcube(cube_t const &bcube_) {bcube = bcube_;}
		unsigned num_roads() const {return roads.size();}
		bool empty() const {return roads.empty();}
		void set_cluster(unsigned id) {cluster_id = id;}
		void register_connected_city(unsigned id) {connected_to.insert(id);}
		set<unsigned> const &get_connected() const {return connected_to;}

		void clear() {
			roads.clear();
			segs.clear();
			plots.clear();
			for (unsigned i = 0; i < 3; ++i) {isecs[i].clear();}
			streetlights.clear();
			parking_lot_mgr.clear();
			tile_blocks.clear();
		}
		bool gen_road_grid(float road_width, float road_spacing) {
			cube_t const &region(bcube); // use our bcube as the region to process
			vector3d const size(region.get_size());
			assert(size.x > 0.0 && size.y > 0.0);
			float const half_width(0.5*road_width), zval(region.z1() + ROAD_HEIGHT);
			float const rx1(region.x1() + half_width), rx2(region.x2() - half_width), ry1(region.y1() + half_width), ry2(region.y2() - half_width); // shrink to include centerlines
			float road_pitch_x(road_width + road_spacing), road_pitch_y(road_pitch_x);
			int const num_x_roads((rx2 - rx1)/road_pitch_x), num_y_roads((ry2 - ry1)/road_pitch_y);
			road_pitch_x = 0.9999*(rx2 - rx1)/num_x_roads; // auto-calculate, round down slightly to avoid FP error
			road_pitch_y = 0.9999*(ry2 - ry1)/num_y_roads;

			// create a grid, for now; crossing roads will overlap
			for (float x = rx1; x < rx2; x += road_pitch_x) {
				roads.emplace_back(point(x, region.y1(), zval), point(x, region.y2(), zval), road_width, true);
			}
			unsigned const num_x(roads.size());

			for (float y = ry1; y < ry2; y += road_pitch_y) {
				roads.emplace_back(point(region.x1(), y, zval), point(region.x2(), y, zval), road_width, false);
			}
			unsigned const num_r(roads.size()), num_y(num_r - num_x);
			if (num_x <= 1 || num_y <= 1) {clear(); return 0;} // not enough space for roads
			bcube.x1() = roads[0      ].x1(); // actual bcube x1 from first x road
			bcube.x2() = roads[num_x-1].x2(); // actual bcube x2 from last  x road
			bcube.y1() = roads[num_x  ].y1(); // actual bcube y1 from first y road
			bcube.y2() = roads[num_r-1].y2(); // actual bcube y2 from last  y road

			// create road segments and intersections
			segs .reserve(num_x*(num_y-1) + (num_x-1)*num_y + 4); // X + Y segments, allocate one extra per side for connectors
			plots.reserve((num_x-1)*(num_y-1));

			if (num_x > 2 && num_y > 2) {
				isecs[0].reserve(4); // 2-way, always exactly 4 at each corner
				isecs[1].reserve(2*((num_x-2) + (num_y-2)) + 4); // 3-way, allocate one extra per side for connectors
				isecs[2].reserve((num_x-2)*(num_y-2)); // 4-way
			}
			for (unsigned x = 0; x < num_x; ++x) {
				for (unsigned y = num_x; y < num_r; ++y) {
					bool const FX(x == 0), FY(y == num_x), LX(x+1 == num_x), LY(y+1 == num_r);
					cube_t const &rx(roads[x]), &ry(roads[y]);
					unsigned const num_conn((!FX) + (!LX) + (!FY) + (!LY));
					if (num_conn < 2) continue; // error?
					uint8_t const conn(((!FX) << 0) | ((!LX) << 1) | ((!FY) << 2) | ((!LY) << 3)); // 1-15
					isecs[num_conn - 2].emplace_back(cube_t(rx.x1(), rx.x2(), ry.y1(), ry.y2(), zval, zval), y, x, conn, false); // intersections
					
					if (!LX) { // skip last y segment
						cube_t const &rxn(roads[x+1]);
						segs.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y1(), ry.y2(), zval, zval), y, false); // y-segments
					}
					if (!LY) { // skip last x segment
						cube_t const &ryn(roads[y+1]);
						segs.emplace_back(cube_t(rx.x1(), rx.x2(), ry.y2(), ryn.y1(), zval, zval), x, true); // x-segments

						if (!LX) { // skip last y segment
							cube_t const &rxn(roads[x+1]);
							plots.push_back(cube_t(rx.x2(), rxn.x1(), ry.y2(), ryn.y1(), zval, zval)); // plots between roads
						}
					}
				} // for y
			} // for x
			return 1;
		}
		void calc_bcube_from_roads() {
			if (roads.empty()) return; // no roads
			bcube = roads.front(); // first road
			for (auto r = roads.begin()+1; r != roads.end(); ++r) {bcube.union_with_cube(*r);} // skip first road
		}
	private:
		int find_conn_int_seg(cube_t const &c, bool dim, bool dir) const {
			float const min_seg_len(1.0*city_params.road_width);

			for (unsigned i = 0; i < segs.size(); ++i) {
				road_seg_t const &s(segs[i]);
				if (s.dim == dim) continue; // not perp dim
				if (s.d[dim][dir] != bcube.d[dim][dir]) continue; // not on edge of road grid
				if (s.d[!dim][1] < c.d[!dim][0] || s.d[!dim][0] > c.d[!dim][1]) continue; // no overlap/projection in other dim
				// c contained in segment in other dim with enough padding (min road length) on each side
				if (c.d[!dim][0] > s.d[!dim][0]+min_seg_len && c.d[!dim][1] < s.d[!dim][1]-min_seg_len) return i; // this is the one we want
				return -1; // partial overlap in other dim, can't split, fail
			} // for i
			return -1; // not found
		}
		int find_conn_intersection(cube_t const &c, bool dim, bool dir, unsigned num_conn) const {
			assert(num_conn == 2 || num_conn == 3); // 2 or 3 way intersections; 4 way can't have more connections
			vector<road_isec_t> const &v(isecs[num_conn-2]);
			for (unsigned i = 0; i < v.size(); ++i) {if (c.intersects_xy(v[i])) return i;}
			return -1; // not found
		}
		void make_4way_int(unsigned int3_ix) { // turn a 3-way intersection into a 4-way intersection for a connector road
			assert(int3_ix < isecs[1].size());
			road_isec_t &isec(isecs[1][int3_ix]); // bbox doesn't change, only conn changes
			isec.make_4way(); // all connected
			isecs[2].push_back(isec); // add as 4-way intersection
			isecs[1][int3_ix] = isecs[1].back(); // remove original 3-way intersection
			isecs[1].pop_back();
		}

		struct road_ixs_t {
			vector<unsigned> seg_ixs, isec_ixs[3][2]; // {2-way, 3-way, 4-way} x {X, Y}
		};
		template<typename T> int search_for_adj(vector<T> const &v, vector<unsigned> const &ixs, cube_t const &bcube, bool dim, bool dir) const {
			for (auto i = ixs.begin(); i != ixs.end(); ++i) {
				assert(*i < v.size());
				cube_t const &c(v[*i]);
				if (c.d[dim][!dir] != bcube.d[dim][dir]) continue; // no shared edge
				if (c.d[!dim][0] != bcube.d[!dim][0] || c.d[!dim][1] != bcube.d[!dim][1]) continue; // no shared edge in other dim
				return *i; // there can be only one
			} // for i
			return -1; // not found
		}
		vector<unsigned> const &get_segs_connecting_to_city(unsigned city) const {
			assert(city < city_to_seg.size());
			return city_to_seg[city];
		}
	public:
		void calc_ix_values(vector<road_network_t> const &road_networks, road_network_t const &global_rn) {
			// now that the segments and intersections are in order, we can fill in the IDs; first, create a mapping from road to segments and intersections
			bool const is_global_rn(&global_rn == this);
			assert(road_to_city.size() == (is_global_rn ? roads.size() : 0));
			vector<road_ixs_t> by_ix(roads.size()); // maps road_ix to list of seg_ix values

			if (is_global_rn) {
				unsigned num_cities(0);
				for (auto r = road_to_city.begin(); r != road_to_city.end(); ++r) {
					for (unsigned d = 0; d < 2; ++d) {if (r->id[d] != CONN_CITY_IX) {max_eq(num_cities, r->id[d]+1);}}
				}
				city_to_seg.resize(num_cities);
			}
			for (unsigned i = 0; i < segs.size(); ++i) {
				unsigned const ix(segs[i].road_ix);
				assert(ix < by_ix.size());
				assert(roads[ix].dim == segs[i].dim);
				by_ix[ix].seg_ixs.push_back(i);
			}
			for (unsigned n = 0; n < 3; ++n) { // {2-way, 3-way, 4-way}
				for (unsigned i = 0; i < isecs[n].size(); ++i) {
					for (unsigned d = 0; d < 2; ++d) { // {x, y}
						int ix(isecs[n][i].rix_xy[d]);

						if (ix < 0) { // global connector road
							ix = decode_neg_ix(ix);
							assert((unsigned)ix < global_rn.roads.size());
							continue; // connector road, nothing else to do here?
						}
						else {
							assert((unsigned)ix < roads.size());
							assert(roads[ix].dim == (d != 0));
							by_ix[ix].isec_ixs[n][d].push_back(i);
						}
					} // for d
				} // for i
			} // for n

			// next, connect segments and intersections together by index, using roads as an acceleration structure
			for (unsigned i = 0; i < segs.size(); ++i) {
				road_seg_t &seg(segs[i]);
				road_ixs_t const &rix(by_ix[seg.road_ix]);

				for (unsigned dir = 0; dir < 2; ++dir) { // dir
					bool found(0);
					int const seg_ix(search_for_adj(segs, rix.seg_ixs, seg, seg.dim, (dir != 0)));
					if (seg_ix >= 0) {assert(seg_ix != i); seg.conn_ix[dir] = seg_ix; seg.conn_type[dir] = TYPE_RSEG; found = 1;} // found segment

					for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
						int const isec_ix(search_for_adj(isecs[n], rix.isec_ixs[n][seg.dim], seg, seg.dim, (dir != 0)));
						if (isec_ix >= 0) {assert(!found); seg.conn_ix[dir] = isec_ix; seg.conn_type[dir] = (TYPE_ISEC2 + n); found = 1;} // found intersection
					}
					if (is_global_rn && !found) { // connection to a city
						seg.conn_type[dir] = TYPE_ISEC3; // always connects to a 3-way intersection within the city
						assert(seg.road_ix < road_to_city.size());
						unsigned const city(road_to_city[seg.road_ix].id[dir]);
						assert(city != CONN_CITY_IX); // internal segments should be connected and not get here
						assert(city < road_networks.size());
						road_network_t const &rn(road_networks[city]);
						vector<unsigned> all_ixs(rn.isecs[1].size());
						for (unsigned n = 0; n < all_ixs.size(); ++n) {all_ixs[n] = n;} // all sequential index values
						int const isec_ix(rn.search_for_adj(rn.isecs[1], all_ixs, seg, seg.dim, (dir != 0)));
						assert(isec_ix >= 0); // must be found
						seg.conn_ix[dir] = isec_ix;
						assert(city < city_to_seg.size());
						city_to_seg[city].push_back(i); // add segment ID
						continue;
					}
					assert(found);
				} // for dir
			} // for i
			for (unsigned n = 0; n < 3; ++n) { // 2-way, 3-way, 4-way
				for (unsigned i = 0; i < isecs[n].size(); ++i) {
					road_isec_t &isec(isecs[n][i]);

					for (unsigned d = 0; d < 4; ++d) { // {-x, +x, -y, +y}
						if (!(isec.conn & (1<<d))) continue; // no connection in this position
						unsigned const dim(d>>1), dir(d&1);
						int const ix(isec.rix_xy[dim]);

						if (ix < 0) { // global connector road
							vector<unsigned> const &seg_ids(global_rn.get_segs_connecting_to_city(city_id));
							assert(!seg_ids.empty());
							int const seg_ix(global_rn.search_for_adj(global_rn.segs, seg_ids, isec, (dim != 0), (dir != 0))); // global conn segment
							assert(seg_ix >= 0); // must be found
							isec.conn_ix[d] = encode_neg_ix(seg_ix);
						}
						else { // local segment
							int const seg_ix(search_for_adj(segs, by_ix[ix].seg_ixs, isec, (dim != 0), (dir != 0))); // always connects to a road segment
							assert(seg_ix >= 0); // must be found
							isec.conn_ix[d] = seg_ix;
						}
					} // for d
				} // for i
			} // for n
		}
		bool check_valid_conn_intersection(cube_t const &c, bool dim, bool dir) const {return (find_conn_int_seg(c, dim, dir) >= 0);}
		void insert_conn_intersection(cube_t const &c, bool dim, bool dir, unsigned grn_rix) { // Note: dim is the dimension of the connector road
			int const seg_id(find_conn_int_seg(c, dim, dir));
			assert(seg_id >= 0 && (unsigned)seg_id < segs.size());
			segs.push_back(segs[seg_id]); // clone the segment first
			road_seg_t &seg(segs[seg_id]);
			assert(seg.road_ix < roads.size() && roads[seg.road_ix].dim != dim); // sanity check
			seg        .d[!dim][1] = c.d[!dim][0]; // low part
			segs.back().d[!dim][0] = c.d[!dim][1]; // high part
			cube_t ibc(seg); // intersection bcube
			ibc.d[!dim][0] = c.d[!dim][0]; // copy width from c
			ibc.d[!dim][1] = c.d[!dim][1];
			uint8_t const conns[4] = {7, 11, 13, 14};
			int const other_rix(encode_neg_ix(grn_rix)); // make negative
			isecs[1].emplace_back(ibc, (dim ? seg.road_ix : (int)other_rix), (dim ? other_rix : (int)seg.road_ix), conns[2*(!dim) + dir], true);
		}
		float create_connector_road(cube_t const &bcube1, cube_t const &bcube2, vector<cube_t> &blockers, road_network_t *rn1, road_network_t *rn2, unsigned city1, unsigned city2,
			heightmap_query_t &hq, float road_width, float conn_pos, bool dim, bool check_only=0)
		{
			bool const dir(bcube1.d[dim][0] < bcube2.d[dim][0]);
			if (dir == 0) {swap(city1, city2);} // make {lo, hi}
			point p1, p2;
			p1.z = bcube1.d[2][1];
			p2.z = bcube2.d[2][1];
			p1[!dim] = p2[!dim] = conn_pos;
			p1[ dim] = bcube1.d[dim][ dir];
			p2[ dim] = bcube2.d[dim][!dir];
			bool const slope((p1.z < p2.z) ^ dir);
			road_t const road(p1, p2, road_width, dim, slope, roads.size());
			float const road_len(road.get_length()), delta_z(road.get_dz()), max_slope(city_params.max_road_slope);
			assert(road_len > 0.0 && delta_z >= 0.0);
			if (delta_z/road_len > max_slope) {assert(check_only); return -1.0;} // slope is too high (split segments will have even higher slopes)
			unsigned const x1(hq.get_x_pos(road.x1())), y1(hq.get_y_pos(road.y1())), x2(hq.get_x_pos(road.x2())), y2(hq.get_y_pos(road.y2()));

			if (check_only) { // only need to do these checks in this case
				// use find_conn_intersection(4)/make_4way_int() to create a new 4-way intersection on one city if one of these fails?
				if (rn1 && !rn1->check_valid_conn_intersection(road, dim,  dir)) return -1.0; // invalid, don't make any changes
				if (rn2 && !rn2->check_valid_conn_intersection(road, dim, !dir)) return -1.0;

				for (auto b = blockers.begin(); b != blockers.end(); ++b) {
					if ((rn1 && b->contains_cube(bcube1)) || (rn2 && b->contains_cube(bcube2))) continue; // skip current cities
					// create an intersection if blocker is a road, and happens to be the same elevation?
					if (b->intersects_xy(road)) return -1.0; // bad intersection, fail
				}
				if (hq.any_underwater(x1, y1, x2+1, y2+1)) return -1.0; // underwater (Note: bounds check is done here)
			}
			if (!check_only) { // create intersections and add blocker
				unsigned const grn_rix(roads.size()); // may be wrong end of connector, but doesn't matter?
				if (rn1) {rn1->insert_conn_intersection(road, dim,  dir, grn_rix);}
				if (rn2) {rn2->insert_conn_intersection(road, dim, !dir, grn_rix);}
				float const blocker_padding(max(city_params.road_spacing, 2.0f*city_params.road_border*max(DX_VAL, DY_VAL))); // use road_spacing?
				blockers.push_back(road);
				blockers.back().expand_by(blocker_padding); // add extra padding
			}
			if (road_len <= city_params.conn_road_seg_len) { // simple single road segment case
				if (!check_only) {
					roads.push_back(road);
					road_to_city.emplace_back(city1, city2);
				}
				return hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope]-ROAD_HEIGHT, road.d[2][!slope]-ROAD_HEIGHT, dim, city_params.road_border, check_only);
			}
			unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
			assert(num_segs > 0 && num_segs < 1000); // sanity check
			float const seg_len(road_len/num_segs);
			assert(seg_len <= city_params.conn_road_seg_len);
			road_t rs(road); // keep d[!dim][0], d[!dim][1], dim, and road_ix
			rs.z1() = road.d[2][slope];
			segments.clear();
			float tot_dz(0.0);

			for (unsigned n = 0; n < num_segs; ++n) {
				rs.d[dim][1] = ((n+1 == num_segs) ? road.d[dim][1] : (rs.d[dim][0] + seg_len)); // make sure it ends exactly at the correct location
				point pos;
				pos[ dim] = rs.d[dim][1];
				pos[!dim] = conn_pos;
				rs.z2()   = hq.get_height_at(pos.x, pos.y) + ROAD_HEIGHT; // terrain height at end of segment
				rs.slope  = (rs.z2() < rs.z1());
				
				if (fabs(rs.get_slope_val()) > max_slope) { // slope is too high, clamp z2 to max allowed value
					if (n+1 == num_segs) return -1.0;
					rs.z2() = rs.z1() + seg_len*max_slope*SIGN(rs.get_dz());
				}
				segments.push_back(rs);
				rs.d[dim][0] = rs.d[dim][1]; rs.z1() = rs.z2(); // shift segment end point
			} // for n
			for (auto s = segments.begin(); s != segments.end(); ++s) {
				if (s->z2() < s->z1()) {swap(s->z2(), s->z1());} // swap zvals if needed
				assert(s->is_normalized());
				tot_dz += hq.flatten_for_road(*s, city_params.road_border, check_only);
				
				if (!check_only) {
					roads.push_back(*s);
					road_to_city.emplace_back(city1, city2); // Note: city index is specified even for internal (non-terminal) roads
				}
			}
			if (!check_only) { // post-flatten pass to fix up dirt at road joints - doesn't help much
				for (auto s = segments.begin(); s != segments.end(); ++s) {hq.flatten_for_road(*s, city_params.road_border, 0, 1);} // decrease_only=1
			}
			return tot_dz; // success
		}
		void create_connector_bend(cube_t const &int_bcube, bool dx, bool dy, unsigned road_ix_x, unsigned road_ix_y) {
			uint8_t const conns[4] = {6, 5, 10, 9};
			isecs[0].emplace_back(int_bcube, road_ix_x, road_ix_y, conns[2*dy + dx], true);
			//blockers.push_back(int_bcube); // ???
		}
		void split_connector_roads(float road_spacing) {
			// Note: here we use segs, maybe 2-way isecs for bends, but not plots
			for (auto r = roads.begin(); r != roads.end(); ++r) {
				bool const d(r->dim), slope(r->slope);
				float const len(r->get_length());
				unsigned const rix(r->road_ix); // not (r - roads.begin())
				if (len <= road_spacing) {segs.emplace_back(*r, rix); continue;} // single segment road
				assert(len > 0.0);
				unsigned const num_segs(ceil(len/road_spacing));
				float const seg_len(len/num_segs), z1(r->d[2][slope]), z2(r->d[2][!slope]); // use fixed-length segments
				assert(seg_len <= road_spacing);
				cube_t c(*r); // start by copying the road's bcube
				
				for (unsigned n = 0; n < num_segs; ++n) {
					c.d[d][1] = ((n+1 == num_segs) ? r->d[d][1] : (c.d[d][0] + seg_len)); // make sure it ends exactly at the correct location
					for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
					if (c.z2() < c.z1()) {swap(c.z2(), c.z1());} // swap zvals if needed
					assert(c.is_normalized());
					segs.emplace_back(c, rix, d, r->slope);
					c.d[d][0] = c.d[d][1]; // shift segment end point
				} // for n
			} // for r
		}
		void gen_tile_blocks() {
			tile_blocks.clear(); // should already be empty?
			tile_to_block_map.clear();
			add_tile_blocks(segs,  tile_to_block_map, TYPE_RSEG);
			add_tile_blocks(plots, tile_to_block_map, TYPE_PLOT);
			for (unsigned i = 0; i < 3; ++i) {add_tile_blocks(isecs[i], tile_to_block_map, (TYPE_ISEC2 + i));}
			//cout << "tile_to_block_map: " << tile_to_block_map.size() << ", tile_blocks: " << tile_blocks.size() << endl;
		}
		void gen_parking_lots(vector<car_t> &cars) {
			parking_lot_mgr.gen_parking(plots, cars, city_id);
			add_tile_blocks(parking_lot_mgr.parks, tile_to_block_map, TYPE_PARK_LOT); // need to do this later, after gen_tile_blocks()
			tile_to_block_map.clear(); // no longer needed
		}
		void add_streetlights() {
			streetlights.clear();
			streetlights.reserve(4*plots.size()); // one on each side of each plot
			float const b(-0.015), a(1.0 - b); // spacing from light pos to plot edge

			for (auto i = plots.begin(); i != plots.end(); ++i) {
				streetlights.emplace_back(point((a*i->x1() + b*i->x2()), (0.75*i->y1() + 0.25*i->y2()), i->z2()), -plus_x); // left   edge one   quarter  up
				streetlights.emplace_back(point((a*i->x2() + b*i->x1()), (0.25*i->y1() + 0.75*i->y2()), i->z2()),  plus_x); // right  edge three quarters up
				streetlights.emplace_back(point((0.25*i->x1() + 0.75*i->x2()), (a*i->y1() + b*i->y2()), i->z2()), -plus_y); // bottom edge three quarters right
				streetlights.emplace_back(point((0.75*i->x1() + 0.25*i->x2()), (a*i->y2() + b*i->y1()), i->z2()),  plus_y); // top    edge one   quarter  right
			}
		}
		void get_road_bcubes(vector<cube_t> &bcubes) const {get_all_bcubes(roads, bcubes);}

		void get_plot_bcubes(vector<cube_t> &bcubes) const { // Note: z-values of cubes indicate building height ranges
			if (plots.empty()) return; // connector road city
			unsigned const start(bcubes.size());
			get_all_bcubes(plots, bcubes);
			vector3d const city_radius(0.5*bcube.get_size());
			point const city_center(bcube.get_cube_center());

			for (auto i = bcubes.begin()+start; i != bcubes.end(); ++i) { // set zvals to control building height range, higher in city center
				point const center(i->get_cube_center());
				float const dx(fabs(center.x - city_center.x)/city_radius.x); // 0 at city center, 1 at city perimeter
				float const dy(fabs(center.y - city_center.y)/city_radius.y); // 0 at city center, 1 at city perimeter
				float const hval(1.0 - max(dx*dx, dy*dy)); // square to give higher weight to larger height ranges
				i->z1() = max(0.0, (hval - 0.25)); // bottom of height range
				i->z2() = min(1.0, (hval + 0.25)); // bottom of height range
			} // for i
		}
		bool check_road_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only=1) const {
			if (roads.empty()) return 0;
			point const query_pos(pos - get_query_xlate());
			if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
			if (check_bcubes_sphere_coll(roads, query_pos, radius, xy_only)) return 1;

			if (include_intersections) { // used for global road network
				for (unsigned i = 0; i < 3; ++i) { // {2-way, 3-way, 4-way}
					if (check_bcubes_sphere_coll(isecs[i], query_pos, radius, xy_only)) return 1;
				}
			}
			return 0;
		}
		void draw(road_draw_state_t &dstate, bool shadow_only) {
			if (empty()) return;
			if (!dstate.check_cube_visible(bcube, 1.0, shadow_only)) return; // VFC/too far

			if (shadow_only) {
				for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
					if (!dstate.check_cube_visible(b->bcube, 1.0, shadow_only)) continue; // VFC/too far
					for (unsigned i = 1; i < 3; ++i) {dstate.draw_stoplights(isecs[i], 1);} // intersections (3-way, 4-way)
				}
			}
			else {
				for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
					if (!dstate.check_cube_visible(b->bcube)) continue; // VFC/too far
					dstate.begin_tile(b->bcube.get_cube_center());
					dstate.draw_road_region(segs,  b->ranges[TYPE_RSEG], b->quads[TYPE_RSEG], TYPE_RSEG); // road segments
					dstate.draw_road_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PLOT], TYPE_PLOT); // plots
					dstate.draw_road_region(parking_lot_mgr.parks, b->ranges[TYPE_PARK_LOT], b->quads[TYPE_PARK_LOT], TYPE_PARK_LOT); // parking lots
					bool const draw_stoplights(dstate.check_cube_visible(b->bcube, 0.16)); // use smaller dist_scale
				
					for (unsigned i = 0; i < 3; ++i) { // intersections (2-way, 3-way, 4-way)
						dstate.draw_road_region(isecs[i], b->ranges[TYPE_ISEC2 + i], b->quads[TYPE_ISEC2 + i], (TYPE_ISEC2 + i));
						if (draw_stoplights && i > 0) {dstate.draw_stoplights(isecs[i], 0);}
					}
				} // for b
			}
			draw_streetlights(dstate, shadow_only);
		}
		void draw_streetlights(road_draw_state_t &dstate, bool shadow_only) const {
			if (streetlights.empty()) return;
			//timer_t t("Draw Streetlights");
			select_texture(WHITE_TEX);
			for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {i->draw(dstate.s, dstate.xlate, shadow_only);}
		}
		void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const { // for now, the only light sources added by the road network are city block streetlights
			for (auto i = streetlights.begin(); i != streetlights.end(); ++i) {i->add_dlight(xlate, lights_bcube);}
		}

		// cars
		static float get_car_lane_offset() {return CAR_LANE_OFFSET*city_params.road_width;}

		bool add_car(car_t &car, rand_gen_t &rgen) const {
			if (segs.empty()) return 0; // no segments to place car on
			vector3d const nom_car_size(city_params.get_car_size());

			for (unsigned n = 0; n < 10; ++n) { // make 10 tries
				unsigned const seg_ix(rgen.rand()%segs.size());
				road_seg_t const &seg(segs[seg_ix]); // chose a random segment
				car.dim   = seg.dim;
				car.dir   = (rgen.rand()&1);
				car.max_speed = rgen.rand_uniform(0.66, 1.0); // add some speed variation
				car.cur_road  = seg.road_ix;
				car.cur_seg   = seg_ix;
				car.cur_road_type = TYPE_RSEG;
				vector3d car_sz(nom_car_size); // {length, width, height} // Note: car models should all be the same size
				car.height = car_sz.z;
				point pos;
				float val1(seg.d[seg.dim][0] + 0.6*car_sz.x), val2(seg.d[seg.dim][1] - 0.6*car_sz.x);
				if (val1 >= val2) continue; // failed, try again (connector road junction?)
				pos[!seg.dim]  = 0.5*(seg.d[!seg.dim][0] + seg.d[!seg.dim][1]); // center of road
				pos[!seg.dim] += ((car.dir ^ car.dim) ? -1.0 : 1.0)*get_car_lane_offset(); // place in right lane
				pos[ seg.dim] = rgen.rand_uniform(val1, val2); // place at random pos in segment
				pos.z = seg.z2() + 0.5*car_sz.z; // place above road surface
				if (seg.dim) {swap(car_sz.x, car_sz.y);}
				car.bcube.set_from_point(pos);
				car.bcube.expand_by(0.5*car_sz);
				assert(get_road_bcube_for_car(car).contains_cube_xy(car.bcube)); // sanity check
				return 1; // success
			} // for n
			return 0; // failed
		}
		void find_car_next_seg(car_t &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
			if (car.cur_road_type == TYPE_RSEG) {
				road_seg_t const &seg(get_car_seg(car));
				car.cur_road_type = seg.conn_type[car.dir];
				car.cur_road      = seg.road_ix;
				car.cur_seg       = seg.conn_ix[car.dir];

				if (!road_to_city.empty()) { // on connector road
					assert(car.cur_road < road_to_city.size());
					unsigned const city_ix(road_to_city[car.cur_road].id[car.dir]);
					
					if (car.in_isect() && city_ix != CONN_CITY_IX) {
						road_network_t const &rn(road_networks[city_ix]);
						vector<road_isec_t> const &isecs(rn.isecs[1]);
						car.cur_city = city_ix;
						assert(car.cur_road_type == TYPE_ISEC3); // must be a 3-way intersection
						assert(car.cur_seg  < isecs.size());
						car.cur_road = isecs[car.cur_seg].rix_xy[!car.dim]; // use the road in the other dim, since it must be within the new city
						assert(car.cur_road < rn.roads.size());
						car.entering_city = 1; // flag so that collision detection works
					}
				}
				assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return; // always within same city, no city update
			}
			road_isec_t const &isec(get_car_isec(car)); // conn_ix: {-x, +x, -y, +y}
			unsigned const orient(car.get_orient());
			int conn_ix(isec.conn_ix[orient]), rix(isec.rix_xy[car.dim]);
			assert(isec.conn & (1<<orient));

			if (conn_ix < 0) { // city connector road case, use global_rn
				assert(rix < 0);
				conn_ix = decode_neg_ix(conn_ix);
				assert((unsigned)conn_ix < global_rn.segs.size());
				rix = global_rn.segs[conn_ix].road_ix;
				car.cur_city = CONN_CITY_IX; // move to global road network
			}
			else {
				assert(rix >= 0); // local road
			}
			car.cur_road = (unsigned)rix;
			car.cur_seg  = (unsigned)conn_ix;
			car.cur_road_type = TYPE_RSEG; // always connects to a road segment
			car.entering_city = 0;
			cube_t const bcube(get_road_bcube_for_car(car, global_rn));

			if (!bcube.intersects_xy(car.bcube)) { // sanity check
				cout << "bad intersection:" << endl << car.str() << endl << "bcube: " << bcube.str() << endl;
				assert(0);
			}
		}
		void update_car(car_t &car, rand_gen_t &rgen, vector<road_network_t> const &road_networks, road_network_t const &global_rn) const {
			assert(car.cur_city == city_id);
			if (car.is_parked()) return; // stopped, no update (for now)

			if (car.stopped_at_light) { // Note: is_isect test is here to allow cars to coast through lights when decel is very low
				bool const was_stopped(car.is_stopped());
				if (!car.in_isect() || get_car_isec(car).can_go_now(car)) {car.stopped_at_light = 0;} // can go now
				else {car.decelerate_fast();}
				if (was_stopped) return; // no update needed
			} else {car.accelerate();}

			cube_t const bcube(get_road_bcube_for_car(car));
			if (!bcube.intersects_xy(car.prev_bcube)) {cout << car.str() << endl << bcube.str() << endl; assert(0);} // sanity check
			bool const dim(car.dim);
			float const road_dz(bcube.get_dz());

			if (road_dz != 0.0) { // car on connector road
				assert(car.cur_road_type == TYPE_RSEG);
				assert(car.cur_city == CONN_CITY_IX);
				bool const slope(get_car_seg(car).slope);
				float const car_pos(0.5*(car.bcube.d[dim][0] + car.bcube.d[dim][1])); // center of car in dim
				float const road_len(bcube.d[dim][1] - bcube.d[dim][0]);
				assert(road_len > TOLERANCE);
				float const t((car_pos - bcube.d[dim][0])/road_len); // car pos along road in (0.0, 1.0)
				float const road_z(bcube.d[2][slope] + t*(bcube.d[2][!slope] - bcube.d[2][slope]));
				float const car_len(car.get_length());
				car.dz = ((slope ^ car.dir) ? 1.0 : -1.0)*road_dz*(car_len/road_len);
				car.bcube.z1() = road_z - 0.5*fabs(car.dz);
				car.bcube.z2() = road_z + 0.5*fabs(car.dz) + car.height;
			}
			else if (car.dz != 0.0) { // car moving from connector road to level city
				float const road_z(bcube.d[2][1]);
				car.dz = 0.0;
				car.bcube.z1() = road_z;
				car.bcube.z2() = road_z + car.height;
			}
			if (car.turn_dir != TURN_NONE) {
				assert(car.in_isect());
				bool const turn_dir(car.turn_dir == TURN_RIGHT); // 0=left, 1=right
				point const car_center(car.get_center()), prev_center(car.prev_bcube.get_cube_center());
				float const car_lane_offset(get_car_lane_offset());
				float const turn_radius((turn_dir ? 0.15 : 0.25)*city_params.road_width); // right turn has smaller radius
				float const isec_center(bcube.get_cube_center()[dim]);
				float const centerline(isec_center + (((car.turn_dir == TURN_LEFT) ^ car.dir) ? -1.0 : 1.0)*car_lane_offset);
				float const prev_val(prev_center[dim]), cur_val(car_center[dim]);
				float const dist_to_turn(fabs(cur_val - centerline));

				if (dist_to_turn < turn_radius) { // turn radius; Note: cars turn around their center points, not their front wheels, which looks odd
					float const dist_from_turn_start(turn_radius - dist_to_turn);
					float const dev(turn_radius - sqrt(turn_radius*turn_radius - dist_from_turn_start*dist_from_turn_start));
					float const new_center(car.turn_val + dev*((turn_dir^car.dir^dim) ? 1.0 : -1.0));
					float const adj(new_center - car_center[!dim]);
					float const frame_dist(p2p_dist_xy(car_center, prev_center)); // total XY distance the car is allowed to move
					car.rot_z = (turn_dir ? -1.0 : 1.0)*(1.0 - CLIP_TO_01(dist_to_turn/turn_radius));
					car.bcube.d[!dim][0] += adj; car.bcube.d[!dim][1] += adj;
					vector3d const move_dir(car.get_center() - prev_center); // total movement from car + turn
					float const move_dist(move_dir.mag());
					
					if (move_dist > TOLERANCE) { // avoid division by zero
						vector3d const delta(move_dir*(frame_dist/move_dist - 1.0)); // overshoot value due to turn
						car.bcube += delta;
					}
				}
				if (min(prev_val, cur_val) <= centerline && max(prev_val, cur_val) > centerline) { // crossed the lane centerline boundary
					car.move_by(centerline - cur_val); // align to lane centerline
					vector3d const car_sz(car.bcube.get_size());
					float const size_adj(0.5*(car_sz[dim] - car_sz[!dim]));
					vector3d expand(zero_vector);
					expand[dim] -= size_adj; expand[!dim] += size_adj;
					car.bcube.expand_by(expand); // fix aspect ratio
					if ((dim == 0) ^ (car.turn_dir == TURN_LEFT)) {car.dir ^= 1;}
					car.dim     ^= 1;
					car.rot_z    = 0.0;
					car.turn_val = 0.0; // reset
					car.turn_dir = TURN_NONE; // turn completed
					car.entering_city = 0;
					road_isec_t const &isec(get_car_isec(car));
					if (isec.conn_ix[car.get_orient()] >= 0) {car.cur_road = isec.rix_xy[car.dim];} // switch to using road_ix in new dim
				}
			}
			if (bcube.contains_cube_xy(car.bcube)) { // in same road seg/int
				assert(get_road_bcube_for_car(car).intersects_xy(car.bcube)); // sanity check
				return; // done
			}
			point car_front(car.get_center());
			car_front[dim] += (car.dir ? 0.375 : -0.375)*car.get_length(); // near the front, so that we can stop at the intersection

			// car crossing the border of this bcube, update state
			if (!bcube.contains_pt_xy_inc_low_edge(car_front)) { // move to another road seg/int
				find_car_next_seg(car, road_networks, global_rn);
				
				if (car.in_isect()) { // moved into an intersection, choose direction
					road_isec_t const &isec(get_car_rn(car, road_networks, global_rn).get_car_isec(car)); // Note: either == city_id, or just moved from global to a new city
					unsigned const orient_in(2*car.dim + (!car.dir)); // invert dir (incoming, not outgoing)
					assert(isec.conn & (1<<orient_in)); // car must come from an enabled orient
					unsigned orients[3]; // {straight, left, right}
					unsigned const conn_left[4] = {3,2,0,1}, conn_right[4] = {2,3,1,0};
					orients[TURN_NONE ] = car.get_orient(); // straight
					orients[TURN_LEFT ] = conn_left [orient_in];
					orients[TURN_RIGHT] = conn_right[orient_in];

					// Note: Could use A* path finding, but it's unlikely to make a visual difference to the casual observer
					while (1) {
						unsigned new_turn_dir(0); // force turn on global conn road 75% of the time to get more cars traveling between cities
						bool const force_turn(isec.is_global_conn_int() && (rgen.rand()&3) != 0);
						int const rval(rgen.rand()%(force_turn ? 2 : 4));
						if      (rval == 0) {new_turn_dir = TURN_LEFT ;} // 25%
						else if (rval == 1) {new_turn_dir = TURN_RIGHT;} // 25%
						else                {new_turn_dir = TURN_NONE ;} // 50%
						if (new_turn_dir == car.front_car_turn_dir && (rgen.rand()%4) != 0) continue; // car in front is too slow, don't turn the same way as it
						if (isec.conn & (1<<orients[new_turn_dir])) {car.turn_dir = new_turn_dir; break;} // success
					} // end while
					assert(isec.conn & (1<<orients[car.turn_dir]));
					car.front_car_turn_dir = TURN_UNSPEC; // reset state now that it's been used
					car.stopped_at_light   = isec.red_or_yellow_light(car);
					if (car.stopped_at_light) {car.decelerate_fast();}
					if (car.turn_dir != TURN_NONE) {car.turn_val = car.get_center()[!dim];} // capture car centerline before the turn
				}
			}
			assert(get_car_rn(car, road_networks, global_rn).get_road_bcube_for_car(car, global_rn).intersects_xy(car.bcube)); // sanity check
		}
		void next_frame() {
			for (unsigned n = 1; n < 3; ++n) { // {2-way, 3-way, 4-way} - Note: 2-way can be skipped
				for (auto i = isecs[n].begin(); i != isecs[n].end(); ++i) {i->next_frame();} // update stoplight state
			}
		}
		static road_network_t const &get_car_rn(car_t const &car, vector<road_network_t> const &road_networks, road_network_t const &global_rn) {
			if (car.cur_city == CONN_CITY_IX) return global_rn;
			assert(car.cur_city < road_networks.size());
			return road_networks[car.cur_city];
		}
		road_seg_t const &get_car_seg(car_t const &car) const {
			assert(car.cur_road_type == TYPE_RSEG);
			assert(car.cur_seg < segs.size());
			return segs[car.cur_seg];
		}
		road_isec_t const &get_car_isec(car_t const &car) const {
			assert(car.in_isect());
			auto const &iv(isecs[car.cur_road_type - TYPE_ISEC2]);
			assert(car.cur_seg < iv.size());
			return iv[car.cur_seg];
		}
		cube_t get_road_bcube_for_car(car_t const &car) const {
			assert(car.cur_road < roads.size()); // Note: generally holds, but not required
			if (car.cur_road_type == TYPE_RSEG) {return get_car_seg(car);}
			return get_car_isec(car);
		}
		cube_t get_road_bcube_for_car(car_t const &car, road_network_t const &global_rn) const { // function variant that works with both global and local roads
			if (car.cur_city == city_id) {return get_road_bcube_for_car(car);}
			else if (car.cur_city == CONN_CITY_IX) {return global_rn.get_road_bcube_for_car(car);}
			else {assert(0); return cube_t();}
		}
	}; // road_network_t

	vector<road_network_t> road_networks; // one per city
	road_network_t global_rn; // connects cities together; no plots
	road_draw_state_t dstate;
	rand_gen_t rgen;

	static float rgen_uniform(float val1, float val2, rand_gen_t &rgen) {return (val1 + (val2 - val1)*rgen.rand_float());}

	void assign_city_clusters() {
		vector<unsigned char> used(road_networks.size(), 0);
		vector<unsigned> pend;
		unsigned cluster_id(0);

		for (unsigned i = 0; i < used.size(); ++i) {
			if (used[i]) continue; // already used
			pend.push_back(i);

			while (!pend.empty()) {
				unsigned const id(pend.back());
				pend.pop_back();
				road_networks[id].set_cluster(cluster_id);
				used[id] = 1;
				set<unsigned> const &conn(road_networks[id].get_connected());

				for (auto c = conn.begin(); c != conn.end(); ++c) {
					assert(*c < used.size());
					if (!used[*c]) {pend.push_back(*c);}
				}
			} // while
			++cluster_id;
		} // for i
		cout << "City clusters: " << cluster_id << endl;
	}

public:
	bool empty() const {return road_networks.empty();}

	cube_t const &get_city_bcube(unsigned city_ix) const {
		if (city_ix == CONN_CITY_IX) {return global_rn.get_bcube();}
		assert(city_ix < road_networks.size());
		return road_networks[city_ix].get_bcube();
	}
	void gen_roads(cube_t const &region, float road_width, float road_spacing) {
		//timer_t timer("Gen Roads"); // ~0.5ms
		road_networks.push_back(road_network_t(region, road_networks.size()));
		if (!road_networks.back().gen_road_grid(road_width, road_spacing)) {road_networks.pop_back(); return;}
		//cout << "Roads: " << road_networks.back().num_roads() << endl;
		road_networks.back().add_streetlights();
	}
	bool connect_two_cities(unsigned city1, unsigned city2, vector<cube_t> &blockers, heightmap_query_t &hq, float road_width) {
		assert(city1 < road_networks.size() && city2 < road_networks.size());
		assert(city1 != city2); // check for self reference
		//cout << "Connect city " << city1 << " and " << city2 << endl;
		road_network_t &rn1(road_networks[city1]), &rn2(road_networks[city2]);
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		assert(!bcube1.intersects_xy(bcube2));
		float const min_edge_dist(4.0*road_width), min_jog(2.0*road_width), half_width(0.5*road_width);
		// Note: cost function should include road length, number of jogs, total elevation change, and max slope

		for (unsigned d = 0; d < 2; ++d) { // try for single segment
			float const shared_min(max(bcube1.d[d][0], bcube2.d[d][0])), shared_max(min(bcube1.d[d][1], bcube2.d[d][1]));
			
			if (shared_max - shared_min > min_edge_dist) { // can connect with single road segment in dim !d, if the terrain in between is passable
				float const val1(shared_min + 0.5*min_edge_dist), val2(shared_max - 0.5*min_edge_dist);
				float best_conn_pos(0.0), best_cost(-1.0);

				for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a straight line
					float const conn_pos(rgen_uniform(val1, val2, rgen)); // chose a random new connection point and try it
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, hq, road_width, conn_pos, !d, 1)); // check_only=1
					if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost;}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Single segment dim: << "d " << cost: " << best_cost << endl;
					global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, city1, city2, hq, road_width, best_conn_pos, !d, 0); // check_only=0; actually make the change
					return 1;
				}
			}
		} // for d
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		bool const dx(center1.x < center2.x), dy(center1.y < center2.y), slope(dx ^ dy);
		cube_t const bc[2] = {bcube1, bcube2};
		
		if ((bc[dx].x1() - bc[!dx].x1() > min_jog) && (bc[dy].y1() - bc[!dy].y1() > min_jog)) {
			// connect with two road segments using a jog: Note: assumes cities are all the same size
			bool const inv_dim(rgen.rand()&1);
			//cout << "Try connect using jog in dim " << inv_dim << endl;
			cube_t bc1_conn(bcube1), bc2_conn(bcube2);
			bc1_conn.d[0][ dx] = bcube2.d[0][!dx];
			bc2_conn.d[0][!dx] = bcube1.d[0][ dx];
			bc1_conn.d[1][ dy] = bcube2.d[1][!dy];
			bc2_conn.d[1][!dy] = bcube1.d[1][ dy];

			for (unsigned d = 0; d < 2; ++d) { // x-then-y vs. y-then-x
				bool const fdim((d != 0) ^ inv_dim); // first segment dim: 0=x first, 1=y first
				bool const range_dir1(fdim ? dy : dx), range_dir2(fdim ? dx : dy);
				cube_t region1(bc1_conn), region2(bc2_conn); // regions of valid jog connectors
				region1.d[ fdim][!range_dir1] = region1.d[ fdim][ range_dir1]; // fixed edge
				region2.d[!fdim][ range_dir2] = region2.d[!fdim][!range_dir2]; // fixed edge
				region1.d[!fdim][0] += min_edge_dist; region1.d[!fdim][1] -= min_edge_dist; // variable edges
				region2.d[ fdim][0] += min_edge_dist; region1.d[ fdim][1] -= min_edge_dist; // variable edges
				float const xmin(region1.d[!fdim][0]), xmax(region1.d[!fdim][1]), ymin(region2.d[fdim][0]), ymax(region2.d[fdim][1]);
				float best_xval(0.0), best_yval(0.0), best_cost(-1.0);
				cube_t best_int_cube;

				for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a single jog
					float xval(rgen_uniform(xmin, xmax, rgen)), yval(rgen_uniform(ymin, ymax, rgen));
					if (!fdim) {swap(xval, yval);}
					float const height(hq.get_height_at(xval, yval) + ROAD_HEIGHT);
					cube_t const int_cube(xval-half_width, xval+half_width, yval-half_width, yval+half_width, height, height); // the candidate intersection point
					bool has_int(0);

					for (auto b = blockers.begin(); b != blockers.end(); ++b) {
						if (b->intersects_xy(int_cube)) {has_int = 1; break;} // bad intersection, fail
					}
					//cout << TXT(dx) << TXT(dy) << TXT(slope) << TXT(fdim) << TXT(range_dir1) << TXT(range_dir2) << TXT(xval) << TXT(yval) << TXT(height) << TXT(has_int) << endl;
					if (has_int) continue; // bad intersection, fail
					float const cost1(global_rn.create_connector_road(bcube1, int_cube, blockers, &rn1, nullptr, city1, CONN_CITY_IX, hq, road_width, (fdim ? xval : yval),  fdim, 1)); // check_only=1
					if (cost1 < 0.0) continue; // bad segment
					if (best_cost > 0.0 && cost1 > best_cost) continue; // bound - early terminate
					float const cost2(global_rn.create_connector_road(int_cube, bcube2, blockers, nullptr, &rn2, CONN_CITY_IX, city2, hq, road_width, (fdim ? yval : xval), !fdim, 1)); // check_only=1
					if (cost2 < 0.0) continue; // bad segment
					float const cost(cost1 + cost2); // Note: cost function will prefer shorter routes
					if (best_cost < 0.0 || cost < best_cost) {best_xval = xval; best_yval = yval; best_int_cube = int_cube; best_cost = cost;}
				} // for n
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Double segment cost: " << best_cost << " " << TXT(best_xval) << TXT(best_yval) << TXT(fdim) << ", int_cube: " << best_int_cube.str() << endl;
					hq.flatten_region_to(best_int_cube, city_params.road_border); // do this first to improve flattening
					unsigned road_ix[2];
					road_ix[ fdim] = global_rn.num_roads();
					global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, city1, CONN_CITY_IX, hq, road_width, (fdim ? best_xval : best_yval),  fdim, 0); // check_only=0
					flatten_op_t const fop(hq.last_flatten_op); // cache for reuse later during decrease_only pass
					road_ix[!fdim] = global_rn.num_roads();
					global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2, CONN_CITY_IX, city2, hq, road_width, (fdim ? best_yval : best_xval), !fdim, 0); // check_only=0
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), road_ix[0], road_ix[1]);
					hq.flatten_sloped_region(fop.x1, fop.y1, fop.x2, fop.y2, fop.z1, fop.z2, fop.dim, fop.border, 0, 1); // decrease_only=1; remove any dirt that the prev road added
					hq.flatten_region_to(best_int_cube, city_params.road_border, 1); // one more pass to fix mesh that was raised above the intersection by a sloped road segment
					return 1;
				}
			} // for d
		}
		return 0;
	}
	void connect_all_cities(float *heightmap, unsigned xsize, unsigned ysize, float road_width, float road_spacing) {
		if (road_width == 0.0 || road_spacing == 0.0) return; // no roads
		unsigned const num_cities(road_networks.size());
		if (num_cities < 2) return; // not cities to connect
		timer_t timer("Connect Cities");
		heightmap_query_t hq(heightmap, xsize, ysize);
		vector<unsigned> is_conn(num_cities, 0); // start with all cities unconnected (0=unconnected, 1=connected, 2=done/connect failed
		unsigned cur_city(0), num_conn(0), num_done(0);
		vector<pair<float, unsigned>> cands;
		vector<cube_t> blockers; // existing cities and connector roads that we want to avoid intersecting

		// gather city blockers
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {
			blockers.push_back(i->get_bcube());
			blockers.back().expand_by(road_spacing); // separate roads by at least this value
		}
		// full cross-product connectivity
		for (unsigned i = 0; i < num_cities; ++i) {
			for (unsigned j = i+1; j < num_cities; ++j) {
				bool const success(connect_two_cities(i, j, blockers, hq, road_width));
				//cout << "Trying to connect city " << i << " to city " << j << ": " << success << endl;
				if (!success) continue;
				road_networks[i].register_connected_city(j);
				road_networks[j].register_connected_city(i);
			} // for j
		} // for i
		assign_city_clusters();
		global_rn.calc_bcube_from_roads();
		global_rn.split_connector_roads(road_spacing);
	}
	void gen_tile_blocks() {
		timer_t timer("Gen Tile Blocks");
		global_rn.gen_tile_blocks(); // must be done first to fill in road_to_city and city_to_seg
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_tile_blocks();}
		global_rn.calc_ix_values(road_networks, global_rn);
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->calc_ix_values(road_networks, global_rn);}
	}
	void gen_parking_lots(vector<car_t> &cars) {
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_parking_lots(cars);}
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {
		global_rn.get_road_bcubes(bcubes); // not sure if this should be included
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_road_bcubes(bcubes);}
	}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const { // Note: no global_rn
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_plot_bcubes(bcubes);}
	}
	bool check_road_sphere_coll(point const &pos, float radius, bool xy_only) const {return global_rn.check_road_sphere_coll(pos, radius, xy_only);}

	void add_city_lights(vector3d const &xlate, cube_t &lights_bcube) const {
		global_rn.add_city_lights(xlate, lights_bcube); // no streetlights, not needed?
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->add_city_lights(xlate, lights_bcube);}
	}
	void draw(int trans_op_mask, vector3d const &xlate, cube_t const &lights_bcube, bool use_dlights, bool shadow_only) { // non-const because qbd is modified
		if (road_networks.empty() && global_rn.empty()) return;

		if (trans_op_mask & 1) { // opaque pass, should be first
			//timer_t timer("Draw Roads");
			fgPushMatrix();
			translate_to(xlate);
			glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
			dstate.pre_draw(xlate, lights_bcube, use_dlights, shadow_only);
			for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate, shadow_only);}
			if (!shadow_only) {global_rn.draw(dstate, shadow_only);} // no stoplights for connector road
			dstate.post_draw();
			glDepthFunc(GL_LESS);
			fgPopMatrix();
		}
		if (trans_op_mask & 2) {dstate.draw_and_clear_light_flares();} // transparent pass; must be done last for alpha blending, and no translate
	}

	// cars
	void next_frame() {
		if (!animate2) return;
		//timer_t timer("Update Stoplights");
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->next_frame();}
		global_rn.next_frame(); // not needed since there are no 3/4-way intersections/stoplights?
	}
	bool add_car(car_t &car, rand_gen_t &rgen) const {
		if (road_networks.empty()) return 0; // no cities to add cars to
		unsigned const city(rgen.rand()%road_networks.size());
		car.cur_city = city;
		return road_networks[city].add_car(car, rgen);
	}
	road_network_t const &get_car_rn(car_t const &car) const {return road_network_t::get_car_rn(car, road_networks, global_rn);}
	void update_car(car_t &car, rand_gen_t &rgen) const {get_car_rn(car).update_car(car, rgen, road_networks, global_rn);}
	cube_t get_road_bcube_for_car(car_t const &car) const {return get_car_rn(car).get_road_bcube_for_car(car);}
	road_isec_t const &get_car_isec(car_t const &car) const {return get_car_rn(car).get_car_isec(car);}
}; // city_road_gen_t


bool car_t::check_collision(car_t &c, city_road_gen_t const &road_gen) {

	if (dir != c.dir) return 0; // traveling on opposite sides of the road
	if (c.dim != dim) return 0; // turning in an intersection, etc.
	float const avg_len(0.5*((bcube.d[dim][1] - bcube.d[dim][0]) + (c.bcube.d[c.dim][1] - c.bcube.d[c.dim][0]))); // average length of the two cars
	float const min_speed(max(0.0f, (min(cur_speed, c.cur_speed) - 0.1f*max_speed))); // relative to max speed of 1.0, clamped to 10% at bottom end for stability
	float const sep_dist(avg_len*(0.25 + 1.11*min_speed)); // 25% to 125% car length, depending on speed
	float const test_dist(0.999*sep_dist); // slightly smaller than separation distance
	cube_t bcube_ext(bcube);
	bcube_ext.d[dim][0] -= test_dist; bcube_ext.d[dim][1] += test_dist; // expand by test_dist distance
	if (!bcube_ext.intersects_xy(c.bcube)) return 0;
	float const front(bcube.d[dim][dir]), c_front(c.bcube.d[c.dim][c.dir]);
	bool const move_c((front < c_front) ^ dir); // move the car that's behind
	// Note: we could slow the car in behind, but that won't work for initial placement collisions when speed == 0
	car_t &cmove(move_c ? c : *this); // the car that will be moved
	car_t const &cstay(move_c ? *this : c); // the car that won't be moved
	//cout << "Collision between " << cmove.str() << " and " << cstay.str() << endl;
	if (cstay.is_stopped()) {cmove.decelerate_fast();} else {cmove.decelerate();}
	float const dist(cstay.bcube.d[dim][!dir] - cmove.bcube.d[dim][dir]); // signed distance between the back of the car in front, and the front of the car in back
	point delta(all_zeros);
	delta[dim] += dist + (cmove.dir ? -sep_dist : sep_dist); // force separation between cars
	cube_t const &bcube(road_gen.get_road_bcube_for_car(cmove));
	if (cstay.max_speed < cmove.max_speed) {cmove.front_car_turn_dir = cstay.turn_dir;} // record the turn dir of this slow car in front of us so we can turn a different way

	if (!bcube.contains_cube_xy(cmove.bcube + delta)) { // moved outside its current road segment bcube
		//if (cmove.bcube == cmove.prev_bcube) {return 1;} // collided, but not safe to move the car (init pos or second collision)
		if (cmove.bcube != cmove.prev_bcube) { // try resetting to last frame's position
			cmove.bcube  = cmove.prev_bcube; // restore prev frame's pos
			return 1; // done
		}
		else { // keep the car from moving outside its current segment (init collision case)
			if (cmove.dir) {max_eq(delta[dim], min(0.0f, 0.999f*(bcube.d[cmove.dim][0] - cmove.bcube.d[cmove.dim][0])));}
			else           {min_eq(delta[dim], max(0.0f, 0.999f*(bcube.d[cmove.dim][1] - cmove.bcube.d[cmove.dim][1])));}
		}
	}
	cmove.bcube += delta;
	return 1;
}


unsigned const NUM_CAR_COLORS = 10;
colorRGBA const car_colors[NUM_CAR_COLORS] = {WHITE, GRAY_BLACK, GRAY, ORANGE, RED, DK_RED, DK_BLUE, DK_GREEN, YELLOW, BROWN};

struct car_model_t {
	string fn;
	int body_mat_id, fixed_color_id;
	float xy_rot, dz, lod_mult; // xy_rot in degrees
	vector<unsigned> shadow_mat_ids;
	car_model_t(string const &fn_, int bmid=-1, int fcid=-1, float rot=0.0, float dz_=0.0, float lm=1.0, vector<unsigned> const &smids=vector<unsigned>()) :
		fn(fn_), body_mat_id(bmid), fixed_color_id(fcid), xy_rot(rot), dz(dz_), lod_mult(lm), shadow_mat_ids(smids) {}
};
unsigned const NUM_CAR_MODELS = 10;
int      const FORCE_MODEL_ID = -1; // -1 disables

car_model_t const car_model_files[NUM_CAR_MODELS] = { // filename, body_material_id, fixed_color_id, xy_rot_angle, delta_z, lod_mult
	car_model_t("../models/cars/sports_car/sportsCar.model3d",        22, -1, 90,  -0.02, 1.0, {20, 22}),
	car_model_t("../models/cars/natla_car/natla_car.obj",             -1,  2, 90,   0.06, 0.5, {1}), // always GRAY
	car_model_t("../models/cars/speedCar/speedCar.obj",               -1,  6, 0,    0.12, 0.5, {4, 5}), // always DK_BLUE
	car_model_t("../models/cars/Lamborghini/Lamborghini.model3d",      2, -1, 180, -0.02, 0.5, {2, 3, 4}),
	car_model_t("../models/cars/GCPD_Police_Car/GCPD_Police_Car.obj", -1,  1, 90,   0.18, 0.2, {0, 1}), // always GRAY_BLACK
	car_model_t("../models/cars/bugatti/bugatti.model3d",              0, -1, 80,  -0.08, 2.0, {0, 4}), // Note: underside disabled for shadows, model is already too many triangles
	car_model_t("../models/cars/Mercedes_Benz/Mercedes-Benz.model3d",  0, -1, 180,  1.00, 0.5, {0, 6, 7}),
	car_model_t("../models/cars/Rio/rio.model3d",                      5, -1, 270,  4.00, 0.5, {1, 5}), // Note: shadow material 1 may be optional
	car_model_t("../models/cars/Soarer/soarer.model3d",                2, -1, 90,   2.00, 0.5, {2, 5, 8}),
	car_model_t("../models/cars/Camaro/camaro2.model3d",              24, -1, 90,   0.10, 0.5, {9, 21, 24}),
	//car_model_t("../models/cars/Bentley/Bentley.model3d",              1, -1, 90,   0.50, 0.5, {1}),
};

class car_manager_t {

	class car_model_loader_t : public model3ds {
		vector<int> models_valid;
		void ensure_models_loaded() {if (empty()) {load_car_models();}}
	public:
		static unsigned num_models() {return NUM_CAR_MODELS;}

		bool is_model_valid(unsigned id) {
			assert(id < NUM_CAR_MODELS);
			ensure_models_loaded(); // I guess we have to load the models here to determine if they're valid
			assert(id < models_valid.size());
			return (models_valid[id] != 0);
		}
		car_model_t const &get_model(unsigned id) const {
			assert(id < NUM_CAR_MODELS);
			return car_model_files[id];
		}
		void load_car_models() {
			models_valid.resize(num_models(), 1); // assume valid

			for (unsigned i = 0; i < num_models(); ++i) {
				string const &fn(get_model(i).fn);
				bool const recalc_normals = 1;

				if (!load_model_file(fn, *this, geom_xform_t(), -1, WHITE, 0, 0.0, recalc_normals, 0, 0, 1)) {
					cerr << "Error: Failed to read model file '" << fn << "'; Skipping this model (will use default box model)." << endl;
					push_back(model3d(fn, tmgr)); // add a placeholder dummy model
					models_valid[i] = 0;
				}
			} // for i
		}
		void draw_car(shader_t &s, vector3d const &pos, cube_t const &car_bcube, vector3d const &dir, colorRGBA const &color, point const &xlate,
			unsigned model_id, bool is_shadow_pass, bool low_detail)
		{
			ensure_models_loaded();
			assert(size() == num_models());
			assert(model_id < size());
			assert(is_model_valid(model_id));
			car_model_t const &model_file(get_model(model_id));
			model3d &model(at(model_id));

			if (!is_shadow_pass && model_file.body_mat_id >= 0) { // use custom color for body material
				material_t &body_mat(model.get_material(model_file.body_mat_id));
				body_mat.ka = body_mat.kd = color;
				//if (model_id == 5) {cout << body_mat.name << endl;}
			}
			model.bind_all_used_tids();
			cube_t const &bcube(model.get_bcube());
			point const orig_camera_pos(camera_pdu.pos);
			camera_pdu.pos += bcube.get_cube_center() - pos - xlate; // required for distance based LOD
			bool const camera_pdu_valid(camera_pdu.valid);
			camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
			// Note: in model space, front-back=z, left-right=x, top-bot=y
			float const sz_scale((car_bcube.get_dx() + car_bcube.get_dy() + car_bcube.get_dz()) / (bcube.get_dx() + bcube.get_dy() + bcube.get_dz()));
			fgPushMatrix();
			translate_to(pos + vector3d(0.0, 0.0, model_file.dz*sz_scale));
			if (fabs(dir.y) > 0.001) {rotate_to_plus_x(dir);}
			else if (dir.x < 0.0) {fgRotate(180.0, 0.0, 0.0, 1.0);}
			fgRotate(TO_DEG*asinf(-dir.z), 0.0, 1.0, 0.0);
			if (model_file.xy_rot != 0.0) {fgRotate(model_file.xy_rot, 0.0, 0.0, 1.0);}
			fgRotate(90.0, 1.0, 0.0, 0.0);
			uniform_scale(sz_scale);
			translate_to(-bcube.get_cube_center()); // cancel out model local translate

			if ((low_detail || is_shadow_pass) && !model_file.shadow_mat_ids.empty()) {
				for (auto i = model_file.shadow_mat_ids.begin(); i != model_file.shadow_mat_ids.end(); ++i) {model.render_material(s, *i, is_shadow_pass);}
			}
			else {
				auto const &unbound_mat(model.get_unbound_material());

				for (unsigned sam_pass = 0; sam_pass < (is_shadow_pass ? 2U : 1U); ++sam_pass) {
					model.render_materials(s, is_shadow_pass, 0, 0, (sam_pass == 1), 3, 3, unbound_mat, rotation_t(),
						nullptr, nullptr, is_shadow_pass, model_file.lod_mult, (is_shadow_pass ? 10.0 : 0.0));
				}
			}
			fgPopMatrix();
			camera_pdu.valid = camera_pdu_valid;
			camera_pdu.pos   = orig_camera_pos;
			select_texture(WHITE_TEX); // reset back to default/untextured
		}
	};

	car_model_loader_t car_model_loader;

	class occlusion_checker_t {
		building_occlusion_state_t state;
	public:
		void set_camera(pos_dir_up const &pdu) {
			if ((display_mode & 0x08) == 0) {state.building_ids.clear(); return;} // testing
			pos_dir_up near_pdu(pdu);
			near_pdu.far_ = 2.0*city_params.road_spacing; // set far clipping plane to one city block
			get_building_occluders(near_pdu, state);
			//cout << "occluders: " << state.building_ids.size() << endl;
		}
		bool is_occluded(cube_t const &c) {
			if (state.building_ids.empty()) return 0;
			float const z(c.z2()); // top edge
			point const corners[4] = {point(c.x1(), c.y1(), z), point(c.x2(), c.y1(), z), point(c.x2(), c.y2(), z), point(c.x1(), c.y2(), z)};
			return check_pts_occluded(corners, 4, state);
		}
	};

	class car_draw_state_t : public draw_state_t {
		quad_batch_draw qbds[3]; // unshadowed, shadowed, AO
		car_model_loader_t &car_model_loader;
		occlusion_checker_t occlusion_checker;
	public:
		car_draw_state_t(car_model_loader_t &car_model_loader_) : car_model_loader(car_model_loader_) {}
		static float get_headlight_dist() {return 3.5*city_params.road_width;} // distance headlights will shine

		colorRGBA get_headlight_color(car_t const &car) const {
			return colorRGBA(1.0, 1.0, (1.0 + 0.8*(fract(1000.0*car.max_speed) - 0.5)), 1.0); // slight yellow-blue tinting using max_speed as a hash
		}
		void pre_draw(vector3d const &xlate_, cube_t const &lights_bcube_, bool use_dlights_, bool shadow_only) {
			//use_bmap = 1; // used only for some car models
			draw_state_t::pre_draw(xlate_, lights_bcube_, use_dlights_, shadow_only, 1); // always_setup_shader=1 (required for model drawing)
			select_texture(WHITE_TEX);
			if (!shadow_only) {occlusion_checker.set_camera(camera_pdu);}
		}
		virtual void draw_unshadowed() {
			qbds[0].draw_and_clear();
			if (qbds[2].empty()) return;
			enable_blend();
			select_texture(BLUR_CENT_TEX);
			qbds[2].draw_and_clear();
			select_texture(WHITE_TEX); // reset back to default/untextured
			disable_blend();
		}
		void add_car_headlights(vector<car_t> const &cars, vector3d const &xlate_, cube_t &lights_bcube) {
			xlate = xlate_; // needed earlier in the flow
			for (auto i = cars.begin(); i != cars.end(); ++i) {add_car_headlights(*i, lights_bcube);}
		}
		void gen_car_pts(car_t const &car, bool include_top, point pb[8], point pt[8]) const {
			point const center(car.get_center());
			cube_t const &c(car.bcube);
			float const z1(center.z - 0.5*car.height), z2(center.z + 0.5*car.height), zmid(center.z + 0.1*car.height), length(car.get_length());
			bool const dim(car.dim), dir(car.dir);
			cube_t top_part(c);
			top_part.d[dim][0] += (dir ? 0.25 : 0.30)*length; // back
			top_part.d[dim][1] -= (dir ? 0.30 : 0.25)*length; // front
			set_cube_pts(c, z1, zmid, dim, dir, pb); // bottom
			if (include_top) {set_cube_pts(top_part, zmid, z2, dim, dir, pt);} // top
			float const sign((dim^dir) ? -1.0 : 1.0);

			if (car.dz != 0.0) { // rotate all points about dim !d
				float const sine_val((dir ? 1.0 : -1.0)*car.dz/length), cos_val(sqrt(1.0 - sine_val*sine_val));
				rotate_pts(center, sine_val, cos_val, dim, 2, pb);
				if (include_top) {rotate_pts(center, sine_val, cos_val, dim, 2, pt);}
			}
			if (car.rot_z != 0.0) { // turning about the z-axis: rot_z of [0.0, 1.0] maps to angles of [0.0, PI/2=90 degrees]
				float const sine_val(sinf(0.5*PI*car.rot_z)), cos_val(sqrt(1.0 - sine_val*sine_val));
				rotate_pts(center, sine_val, cos_val, 0, 1, pb);
				if (include_top) {rotate_pts(center, sine_val, cos_val, 0, 1, pt);}
			}
		}
		void draw_car(car_t const &car, bool shadow_only, bool is_dlight_shadows) { // Note: all quads
			if (is_dlight_shadows) {
				if (!dist_less_than(camera_pdu.pos, car.get_center(), 0.6*camera_pdu.far_)) return;
				cube_t bcube(car.bcube);
				bcube.expand_by(0.1*car.height);
				if (bcube.contains_pt(camera_pdu.pos)) return; // don't self-shadow
			}
			if (!check_cube_visible(car.bcube, (shadow_only ? 0.0 : 0.75))) return; // dist_scale=0.75
			point const center(car.get_center());
			begin_tile(center); // enable shadows
			assert(car.color_id < NUM_CAR_COLORS);
			colorRGBA const &color(car_colors[car.color_id]);
			float const dist_val(p2p_dist(camera_pdu.pos, (center + xlate))/get_draw_tile_dist());
			bool const draw_top(dist_val < 0.25), dim(car.dim), dir(car.dir);
			float const sign((dim^dir) ? -1.0 : 1.0);
			point pb[8], pt[8]; // bottom and top sections
			gen_car_pts(car, draw_top, pb, pt);

			if ((shadow_only || dist_val < 0.05) && car_model_loader.is_model_valid(car.model_id)) {
				if (!shadow_only && occlusion_checker.is_occluded(car.bcube + xlate)) return; // only check occlusion for expensive car models
				vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
				car_model_loader.draw_car(s, center, car.bcube, front_n, color, xlate, car.model_id, shadow_only, (dist_val > 0.035));
			}
			else { // draw simple 1-2 cube model
				quad_batch_draw &qbd(qbds[emit_now]);
				color_wrapper cw; cw.set_c4(color);
				draw_cube(qbd, dim, dir, cw, center, pb); // bottom
				if (draw_top) {draw_cube(qbd, dim, dir, cw, center, pt);} // top
				if (emit_now) {qbds[1].draw_and_clear();} // shadowed (only emit when tile changes?)
			}
			if (dist_val < 0.04 && fabs(car.dz) < 0.01) { // add AO planes when close to the camera and on a level road
				float const length(car.get_length());
				point pao[4];
				
				for (unsigned i = 0; i < 4; ++i) {
					point &v(pao[i]);
					v = pb[i] - center;
					v[ dim] += 0.1*length*SIGN(v[ dim]); // increase length slightly
					v[!dim] += 0.1*length*SIGN(v[!dim]); // increase width  slightly
					v   += center;
					v.z += 0.02*car.height; // shift up slightly to avoid z-fighting
				}
				/*if (!is_night()) { // daytime, adjust shadow to match sun pos
					vector3d const sun_dir(0.5*length*(center - get_sun_pos()).get_norm());
					vector3d const offset(sun_dir.x, sun_dir.y, 0.0);
					for (unsigned i = 0; i < 4; ++i) {pao[i] += offset;} // problems: double shadows, non-flat surfaces, buildings, texture coords/back in center, non-rectangular
				}*/
				qbds[2].add_quad_pts(pao, colorRGBA(0, 0, 0, 0.9), plus_z);
			}
			if (dist_val > 0.3)  return; // to far - no lights to draw
			if (shadow_only || car.is_parked()) return; // no lights when parked, or in shadow pass
			vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
			unsigned const lr_xor(((camera_pdu.pos[!dim] - xlate[!dim]) - center[!dim]) < 0.0);

			if (is_night() && dist_val < 0.3) { // night time headlights
				colorRGBA const hl_color(get_headlight_color(car));

				for (unsigned d = 0; d < 2; ++d) { // L, R
					unsigned const lr(d ^ lr_xor ^ 1);
					point const pos((lr ? 0.2 : 0.8)*(0.2*pb[0] + 0.8*pb[4]) + (lr ? 0.8 : 0.2)*(0.2*pb[1] + 0.8*pb[5]));
					add_light_flare(pos, front_n, hl_color, 2.0, 0.65*car.height); // pb 0,1,4,5
				}
			}
			if ((car.is_almost_stopped() || car.stopped_at_light) && dist_val < 0.2) { // brake lights
				for (unsigned d = 0; d < 2; ++d) { // L, R
					unsigned const lr(d ^ lr_xor);
					point const pos((lr ? 0.2 : 0.8)*(0.2*pb[2] + 0.8*pb[6]) + (lr ? 0.8 : 0.2)*(0.2*pb[3] + 0.8*pb[7]));
					add_light_flare(pos, -front_n, RED, 1.0, 0.5*car.height); // pb 2,3,6,7
				}
			}
			if (car.turn_dir != TURN_NONE && dist_val < 0.1) { // turn signals
				float const ts_period = 1.5; // in seconds
				double const time(fract((tfticks + 1000.0*car.max_speed)/(ts_period*TICKS_PER_SECOND))); // use car max_speed as seed to offset time base
				
				if (time > 0.5) { // flash on and off
					bool const tdir((car.turn_dir == TURN_LEFT) ^ dim ^ dir); // R=1,2,5,6 or L=0,3,4,7
					vector3d const side_n(cross_product((pb[6] - pb[2]), (pb[1] - pb[2])).get_norm()*sign*(tdir ? 1.0 : -1.0));

					for (unsigned d = 0; d < 2; ++d) { // B, F
						point const pos(0.3*pb[tdir ? (d ? 1 : 2) : (d ? 0 : 3)] + 0.7*pb[tdir ? (d ? 5 : 6) : (d ? 4 : 7)]);
						add_light_flare(pos, (side_n + (d ? 1.0 : -1.0)*front_n).get_norm(), colorRGBA(1.0, 0.75, 0.0, 1.0), 1.5, 0.3*car.height); // normal points out 45 degrees
					}
				}
			}
		}
		void add_car_headlights(car_t const &car, cube_t &lights_bcube) {
			if (car.is_parked()) return; // no lights when parked
			float const headlight_dist(get_headlight_dist());
			cube_t bcube(car.bcube);
			bcube.expand_by(headlight_dist);
			if (!lights_bcube.contains_cube_xy(bcube))   return; // not contained within the light volume
			if (!camera_pdu.cube_visible(bcube + xlate)) return; // VFC
			float const sign((car.dim^car.dir) ? -1.0 : 1.0);
			point pb[8], pt[8]; // bottom and top sections
			gen_car_pts(car, 0, pb, pt); // draw_top=0
			vector3d const front_n(cross_product((pb[5] - pb[1]), (pb[0] - pb[1])).get_norm()*sign);
			vector3d const dir((0.5*front_n - 0.5*plus_z).get_norm()); // point slightly down
			colorRGBA const color(get_headlight_color(car));
			float const beamwidth = 0.08;
			min_eq(lights_bcube.z1(), bcube.z1());
			max_eq(lights_bcube.z2(), bcube.z2());

			for (unsigned d = 0; d < 2; ++d) { // L, R
				point const pos((d ? 0.2 : 0.8)*(0.2*pb[0] + 0.8*pb[4]) + (d ? 0.8 : 0.2)*(0.2*pb[1] + 0.8*pb[5]));
				dl_sources.push_back(light_source(headlight_dist, pos, pos, color, 1, dir, beamwidth));
			}
		}
	}; // car_draw_state_t

	struct car_block_t {
		unsigned start, cur_city;
		car_block_t(unsigned s, unsigned c) : start(s), cur_city(c) {}
	};

	city_road_gen_t const &road_gen;
	vector<car_t> cars;
	vector<car_block_t> car_blocks;
	car_draw_state_t dstate;
	rand_gen_t rgen;
	vector<unsigned> entering_city;

public:
	car_manager_t(city_road_gen_t const &road_gen_) : road_gen(road_gen_), dstate(car_model_loader) {}
	bool empty() const {return cars.empty();}
	
	void clear() {
		cars.clear();
		car_blocks.clear();
	}
	void init_cars(unsigned num) {
		if (num == 0) return;
		timer_t timer("Init Cars");
		cars.reserve(num);
		
		for (unsigned n = 0; n < num; ++n) {
			car_t car;
			if (!road_gen.add_car(car, rgen)) continue;
			cars.push_back(car);
		} // for n
		cout << "Dynamic Cars: " << cars.size() << endl;
	}
	void add_parked_cars(vector<car_t> const &new_cars) {
		cars.insert(cars.end(), new_cars.begin(), new_cars.end());
	}
	void finalize_cars() {
		unsigned const num_models(car_model_loader.num_models());

		for (auto i = cars.begin(); i != cars.end(); ++i) {
			int fixed_color(-1);

			if (num_models > 0) {
				if (FORCE_MODEL_ID >= 0) {i->model_id = FORCE_MODEL_ID;}
				else {i->model_id = ((num_models > 1) ? (rgen.rand() % num_models) : 0);}
				fixed_color = car_model_loader.get_model(i->model_id).fixed_color_id;
			}
			i->color_id = ((fixed_color >= 0) ? fixed_color : (rgen.rand() % NUM_CAR_COLORS));
			assert(i->is_valid());
		} // for i
		cout << "Total Cars: " << cars.size() << endl;
	}
	void next_frame(float car_speed) {
		if (cars.empty() || !animate2) return;
		//timer_t timer("Update Cars"); // 10K cars = 1.7ms / 2K cars = 0.2ms
		sort(cars.begin(), cars.end(), comp_car_road_then_pos(dstate.xlate)); // sort by city/road/position for intersection tests and tile shadow map binds
		entering_city.clear();
		car_blocks.clear();
		float const speed(0.001*car_speed*fticks);
		//unsigned num_on_conn_road(0);
		
		for (auto i = cars.begin(); i != cars.end(); ++i) { // move cars
			if (car_blocks.empty() || i->cur_city != car_blocks.back().cur_city) {car_blocks.emplace_back((i - cars.begin()), i->cur_city);}
			if (i->is_parked()) continue; // no update for parked cars
			i->move(speed);
			if (i->entering_city) {entering_city.push_back(i - cars.begin());} // record for use in collision detection
			if (!i->stopped_at_light && i->in_isect()) {road_gen.get_car_isec(*i).stoplight.mark_blocked(i->dim, i->dir);} // blocking intersection
		}
		car_blocks.emplace_back(cars.size(), 0); // add terminator

		for (auto i = cars.begin(); i != cars.end(); ++i) { // collision detection
			if (i->is_parked()) continue; // no collisions for parked cars
			bool const on_conn_road(i->cur_city == CONN_CITY_IX);

			for (auto j = i+1; j != cars.end(); ++j) { // check for collisions with cars on the same road (can't test seg because they can be on diff segs but still collide)
				if (i->cur_city != j->cur_city || i->cur_road != j->cur_road) break; // different cities or roads
				if (!on_conn_road && i->cur_road_type == j->cur_road_type && abs((int)i->cur_seg - (int)j->cur_seg) > (on_conn_road ? 1 : 0)) break; // diff road segs or diff isects
				i->check_collision(*j, road_gen);
			}
			if (on_conn_road) { // on connector road, check before entering intersection to a city
				for (auto ix = entering_city.begin(); ix != entering_city.end(); ++ix) {
					if (*ix != (i - cars.begin())) {i->check_collision(cars[*ix], road_gen);}
				}
				//++num_on_conn_road;
			}
		} // for i
		for (auto i = cars.begin(); i != cars.end(); ++i) {road_gen.update_car(*i, rgen);} // run update logic
		//cout << TXT(cars.size()) << TXT(entering_city.size()) << TXT(in_isects.size()) << TXT(num_on_conn_road) << endl; // TESTING
	}
	void draw(int trans_op_mask, vector3d const &xlate, cube_t const &lights_bcube, bool use_dlights, bool shadow_only) {
		if (cars.empty()) return;

		if (trans_op_mask & 1) { // opaque pass, should be first
			bool const is_dlight_shadows(shadow_only && xlate == zero_vector); // not the best way to test for this, should make shadow_only 3-valued
			if (is_dlight_shadows && !city_params.car_shadows) return;
			bool const only_parked(shadow_only && !is_dlight_shadows); // sun/moon shadows are precomputed and cached, so only include static objects such as parked cars
			//timer_t timer(string("Draw Cars") + (shadow_only ? " Shadow" : "")); // 10K cars = 1.5ms / 2K cars = 0.33ms
			dstate.xlate = xlate;
			fgPushMatrix();
			translate_to(xlate);
			dstate.pre_draw(xlate, lights_bcube, use_dlights, shadow_only);

			for (auto cb = car_blocks.begin(); cb+1 < car_blocks.end(); ++cb) {
				if (!camera_pdu.cube_visible(road_gen.get_city_bcube(cb->cur_city) + xlate)) continue; // city not visible - skip
				unsigned const end((cb+1)->start);
				assert(end <= cars.size());
				
				for (unsigned c = cb->start; c != end; ++c) {
					if (only_parked && !cars[c].is_parked()) continue; // skip non-parked cars
					dstate.draw_car(cars[c], shadow_only, is_dlight_shadows);
				}
			} // for cb
			dstate.post_draw();
			fgPopMatrix();
		}
		if ((trans_op_mask & 2) && !shadow_only) {dstate.draw_and_clear_light_flares();} // transparent pass; must be done last for alpha blending, and no translate
	}
	void add_car_headlights(vector3d const &xlate, cube_t &lights_bcube) {dstate.add_car_headlights(cars, xlate, lights_bcube);}
	void free_context() {car_model_loader.free_context();}
}; // car_manager_t


struct cmp_light_source_sz_dist {
	point const &cpos;
	cmp_light_source_sz_dist(point const &p) : cpos(p) {}
	float get_value(light_source const &s) const {return s.get_beamwidth()*s.get_radius()*s.get_radius()/p2p_dist_sq(s.get_pos(), cpos);}
	bool operator()(light_source const &a, light_source const &b) const {return (get_value(a) > get_value(b));} // sort largest/closest to smallest/furthest
};

void sort_lights_by_dist_size(vector<light_source> &lights, point const &cpos) {
	stable_sort(lights.begin(), lights.end(), cmp_light_source_sz_dist(cpos));
}

void filter_dlights_to(vector<light_source> &lights, unsigned max_num, point const &cpos) {
	if (lights.size() <= max_num) return;
	sort_lights_by_dist_size(lights, cpos);
	lights.resize(max_num); // remove lowest scoring lights
}

struct city_smap_manager_t {
	void setup_shadow_maps(vector<light_source> &light_sources, point const &cpos) {
		unsigned const num_smaps(min(light_sources.size(), min(city_params.max_shadow_maps, MAX_DLIGHT_SMAPS)));
		dl_smap_enabled = 0;
		if (!enable_dlight_shadows || shadow_map_sz == 0 || num_smaps == 0) return;
		sort_lights_by_dist_size(light_sources, cpos);
		unsigned num_used(0);
		
		// Note: slow to recreate shadow maps every frame, but most lights are either dynamic (headlights) or include dynamic shadow casters (cars) and need to be updated every frame anyway
		for (auto i = light_sources.begin(); i != light_sources.end() && num_used < num_smaps; ++i) {
			if (!i->is_very_directional()) continue; // not a spotlight
			//if (!city_params.car_shadows && i->is_dynamic()) continue; // skip headlights (optimization)
			dl_smap_enabled |= i->setup_shadow_map(CITY_LIGHT_FALLOFF);
			++num_used;
		} // for i
	}
	void clear_all_smaps(vector<light_source> &light_sources) {
		for (auto i = light_sources.begin(); i != light_sources.end(); ++i) {i->release_smap();}
	}
};

class city_gen_t : public city_plot_gen_t {

	city_road_gen_t road_gen;
	car_manager_t car_manager;
	city_smap_manager_t city_smap_manager;
	cube_t lights_bcube;
	float light_radius_scale;

public:
	city_gen_t() : car_manager(road_gen), lights_bcube(all_zeros), light_radius_scale(1.0) {}

	bool gen_city(city_params_t const &params, cube_t &cities_bcube) {
		unsigned x1(0), y1(0), x2(0), y2(0);
		if (!find_best_city_location(params.city_size_min, params.city_size_min, params.city_size_max, params.city_size_max,
			params.city_border, params.slope_width, params.num_samples, x1, y1, x2, y2)) return 0;
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		if (params.roads_enabled()) {road_gen.gen_roads(pos_range, params.road_width, params.road_spacing);}
		return 1;
	}
	void gen_cities(city_params_t const &params) {
		if (params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		{ // open a scope
			timer_t t("Choose City Location");
			for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params, cities_bcube);}
		}
		bool const is_const_zval(cities_bcube.z1() == cities_bcube.z2());
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube, is_const_zval);}
		road_gen.connect_all_cities(heightmap, xsize, ysize, params.road_width, params.road_spacing);
		road_gen.gen_tile_blocks();
		car_manager.init_cars(city_params.num_cars);
	}
	void gen_details() {
		if (road_gen.empty() || car_manager.empty()) return; // nothing to do - no roads or cars
		// generate parking lots
		vector<car_t> parked_cars;
		road_gen.gen_parking_lots(parked_cars);
		car_manager.add_parked_cars(parked_cars);
		car_manager.finalize_cars();
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_road_bcubes(bcubes);}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_plot_bcubes(bcubes);}

	bool check_city_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (check_plot_sphere_coll(pos, radius, xy_only)) return 1;
		return road_gen.check_road_sphere_coll(pos, radius, xy_only);
	}
	void next_frame() {
		road_gen.next_frame(); // update stoplights
		car_manager.next_frame(city_params.car_speed);
	}
	void draw(bool shadow_only, int reflection_pass, int trans_op_mask, vector3d const &xlate) { // for now, there are only roads
		bool const use_dlights(is_night());
		if (reflection_pass == 0) {road_gen.draw(trans_op_mask, xlate, lights_bcube, use_dlights, shadow_only);} // roads don't cast shadows and aren't reflected in water, but stoplights cast shadows
		car_manager.draw(trans_op_mask, xlate, lights_bcube, use_dlights, shadow_only);
		// Note: buildings are drawn through draw_buildings()
	}
	cube_t setup_city_lights(vector3d const &xlate) {
		bool const prev_had_lights(!dl_sources.empty());
		city_smap_manager.clear_all_smaps(dl_sources);
		clear_dynamic_lights();
		lights_bcube.set_to_zeros();
		if (!is_night() && !prev_had_lights) return lights_bcube; // only have lights at night
		float const light_radius(1.0*light_radius_scale*get_tile_smap_dist()); // distance from the camera where headlights and streetlights are drawn
		point const cpos(camera_pdu.pos - xlate);
		lights_bcube = cube_t(cpos);
		lights_bcube.expand_by(light_radius);
		lights_bcube.z1() =  FLT_MAX;
		lights_bcube.z2() = -FLT_MAX;
		car_manager.add_car_headlights(xlate, lights_bcube);
		road_gen.add_city_lights(xlate, lights_bcube);
		//cout << "dlights: " << dl_sources.size() << ", bcube: " << lights_bcube.str() << endl;
		unsigned const max_dlights(min(1024U, city_params.max_lights)); // Note: should be <= the value used in upload_dlights_textures()

		if (dl_sources.size() > max_dlights) {
			if (dl_sources.size() > 4*max_dlights) { // too many lights, reduce light radius for next frame
				light_radius_scale *= 0.95;
				cout << "Too many city lights: " << dl_sources.size() << ". Reducing light_radius_scale to " << light_radius_scale << endl;
			}
			filter_dlights_to(dl_sources, max_dlights, cpos);
		}
		city_smap_manager.setup_shadow_maps(dl_sources, cpos);
		add_dynamic_lights_city(lights_bcube);
		upload_dlights_textures(lights_bcube);
		return lights_bcube;
	}
	void free_context() {car_manager.free_context();}
}; // city_gen_t

city_gen_t city_gen;


bool parse_city_option(FILE *fp) {return city_params.read_option(fp);}
bool have_cities() {return city_params.enabled();}
float get_road_max_len() {return city_params.road_spacing;}

void gen_cities(float *heightmap, unsigned xsize, unsigned ysize) {
	if (!have_cities()) return; // nothing to do
	city_gen.init(heightmap, xsize, ysize); // only need to call once for any given heightmap
	city_gen.gen_cities(city_params);
}
void gen_city_details() {city_gen.gen_details();} // called after gen_buildings()
void get_city_road_bcubes(vector<cube_t> &bcubes) {city_gen.get_all_road_bcubes(bcubes);}
void get_city_plot_bcubes(vector<cube_t> &bcubes) {city_gen.get_all_plot_bcubes(bcubes);}
void next_city_frame() {city_gen.next_frame();}
void draw_cities(bool shadow_only, int reflection_pass, int trans_op_mask, vector3d const &xlate) {city_gen.draw(shadow_only, reflection_pass, trans_op_mask, xlate);}
cube_t setup_city_lights(vector3d const &xlate) {return city_gen.setup_city_lights(xlate);}

bool check_city_sphere_coll(point const &pos, float radius) {
	if (!have_cities()) return 0;
	point center(pos);
	if (world_mode == WMODE_INF_TERRAIN) {center += vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0);} // apply xlate for all static objects
	return city_gen.check_city_sphere_coll(center, radius);
}
bool check_valid_scenery_pos(point const &pos, float radius) {
	if (check_buildings_sphere_coll(pos, radius, 1, 1)) return 0; // apply_tt_xlate=1, xy_only=1
	if (check_city_sphere_coll(pos, radius)) return 0;
	return 1;
}
void free_city_context() {city_gen.free_context();}

