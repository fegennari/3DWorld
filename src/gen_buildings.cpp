// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "file_utils.h"
#include "buildings.h"
#include "mesh.h"

using std::string;

unsigned const MAX_CYLIN_SIDES = 36;
float const WIND_LIGHT_ON_RAND = 0.08;

extern bool start_in_inf_terrain, no_store_model_textures_in_memory;
extern int rand_gen_index, display_mode;
extern point sun_pos;

void get_all_model_bcubes(vector<cube_t> &bcubes); // from model3d.h


struct tid_nm_pair_t {

	int tid, nm_tid; // Note: assumes each tid has only one nm_tid
	float tscale_x, tscale_y, txoff, tyoff;

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale_x(1.0), tscale_y(1.0), txoff(0.0), tyoff(0.0) {}
	tid_nm_pair_t(int tid_, int nm_tid_, float tx, float ty, float xo=0.0, float yo=0.0) : tid(tid_), nm_tid(nm_tid_), tscale_x(tx), tscale_y(ty), txoff(xo), tyoff(yo) {}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}
	bool operator==(tid_nm_pair_t const &t) const {return (tid == t.tid && nm_tid == t.nm_tid && tscale_x == t.tscale_x && tscale_y == t.tscale_y);}
	colorRGBA get_avg_color() const {return texture_color(tid);}
	tid_nm_pair_t get_scaled_version(float scale) const {return tid_nm_pair_t(tid, nm_tid, scale*tscale_x, scale*tscale_y);}

	void set_gl() const {
		select_texture(tid);
		select_multitex(((nm_tid < 0) ? FLAT_NMAP_TEX : nm_tid), 5);
	}
};

struct building_tex_params_t {
	tid_nm_pair_t side_tex, roof_tex;
};

struct color_range_t {

	float grayscale_rand;
	colorRGBA cmin, cmax; // alpha is unused?
	color_range_t() : grayscale_rand(0.0), cmin(WHITE), cmax(WHITE) {}

	void gen_color(colorRGBA &color, rand_gen_t &rgen) const {
		if (cmin == cmax) {color = cmin;} // single exact color
		else {UNROLL_4X(color[i_] = rgen.rand_uniform(cmin[i_], cmax[i_]);)}
		if (grayscale_rand > 0.0) {
			float const v(grayscale_rand*rgen.rand_float());
			UNROLL_3X(color[i_] += v;)
		}
	}
};

struct building_mat_t : public building_tex_params_t {

	bool no_city, add_windows, add_wind_lights;
	unsigned min_levels, max_levels, min_sides, max_sides;
	float place_radius, max_delta_z, max_rot_angle, min_level_height, min_alt, max_alt, house_prob, house_scale_min, house_scale_max;
	float split_prob, cube_prob, round_prob, asf_prob, min_fsa, max_fsa, min_asf, max_asf, wind_xscale, wind_yscale, wind_xoff, wind_yoff;
	cube_t pos_range, prev_pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color;
	colorRGBA window_color;

	building_mat_t() : no_city(0), add_windows(0), add_wind_lights(0), min_levels(1), max_levels(1), min_sides(4), max_sides(4), place_radius(0.0),
		max_delta_z(0.0), max_rot_angle(0.0), min_level_height(0.0), min_alt(-1000), max_alt(1000), house_prob(0.0), house_scale_min(1.0), house_scale_max(1.0),
		split_prob(0.0), cube_prob(1.0), round_prob(0.0), asf_prob(0.0), min_fsa(0.0), max_fsa(0.0), min_asf(0.0), max_asf(0.0), wind_xscale(1.0),
		wind_yscale(1.0), wind_xoff(0.0), wind_yoff(0.0), pos_range(-100,100,-100,100,0,0), prev_pos_range(all_zeros), sz_range(1,1,1,1,1,1), window_color(GRAY) {}
	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0);}
	float gen_size_scale(rand_gen_t &rgen) const {return ((house_scale_min == house_scale_max) ? house_scale_min : rgen.rand_uniform(house_scale_min, house_scale_max));}

	void update_range(vector3d const &range_translate) {
		if (place_radius > 0.0) { // clip range to place_radius
			point const center(pos_range.get_cube_center());
			
			for (unsigned d = 0; d < 2; ++d) { // x,y
				max_eq(pos_range.d[d][0], (center[d] - place_radius));
				min_eq(pos_range.d[d][1], (center[d] + place_radius));
			}
		}
		pos_range += range_translate;
	}
	void set_pos_range(cube_t const &new_pos_range) {
		prev_pos_range = pos_range;
		pos_range = new_pos_range;
	}
	void restore_prev_pos_range() {
		if (!prev_pos_range.is_all_zeros()) {pos_range = prev_pos_range;}
	}
};

struct building_params_t {

	bool flatten_mesh, has_normal_map, tex_mirror, tex_inv_y, tt_only, infinite_buildings;
	unsigned num_place, num_tries, cur_prob;
	float ao_factor, sec_extra_spacing;
	float window_width, window_height, window_xspace, window_yspace; // windows
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix, mat_gen_ix_city, mat_gen_ix_nocity; // {any, city_only, non_city}

	building_params_t(unsigned num=0) : flatten_mesh(0), has_normal_map(0), tex_mirror(0), tex_inv_y(0), tt_only(0), num_place(num), num_tries(10),
		cur_prob(1), ao_factor(0.0), sec_extra_spacing(0.0), window_width(0.0), window_height(0.0), window_xspace(0.0), window_yspace(0.0), range_translate(zero_vector) {}
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
	bool windows_enabled  () const {return (window_width > 0.0 && window_height > 0.0 && window_xspace > 0.0 && window_yspace);} // all must be specified as nonzero
	bool gen_inf_buildings() const {return (infinite_buildings && world_mode == WMODE_INF_TERRAIN);} // TODO: for future use
	float get_window_width_fract () const {assert(windows_enabled()); return window_width /(window_width  + window_xspace);}
	float get_window_height_fract() const {assert(windows_enabled()); return window_height/(window_height + window_yspace);}
	float get_window_tx() const {assert(windows_enabled()); return 1.0/(window_width  + window_xspace);}
	float get_window_ty() const {assert(windows_enabled()); return 1.0/(window_height + window_yspace);}
	
	void add_cur_mat() {
		unsigned const mat_ix(materials.size());
		
		for (unsigned n = 0; n < cur_prob; ++n) { // add more references to this mat for higher probability
			mat_gen_ix.push_back(mat_ix);
			(cur_mat.no_city ? mat_gen_ix_nocity : mat_gen_ix_city).push_back(mat_ix);
		}
		materials.push_back(cur_mat);
		materials.back().update_range(range_translate);
		has_normal_map |= cur_mat.has_normal_map();
	}
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
	unsigned choose_rand_mat(rand_gen_t &rgen, bool city_only, bool non_city_only) const {
		vector<unsigned> const &mat_ix_list(get_mat_list(city_only, non_city_only));
		assert(!mat_ix_list.empty());
		return mat_ix_list[rgen.rand()%mat_ix_list.size()];
	}
	void set_pos_range(cube_t const &pos_range) {
		//cout << "pos_range: " << pos_range.str() << endl;
		cur_mat.set_pos_range(pos_range);
		for (auto i = materials.begin(); i != materials.end(); ++i) {i->set_pos_range(pos_range);}
	}
	void restore_prev_pos_range() {
		cur_mat.restore_prev_pos_range();
		for (auto i = materials.begin(); i != materials.end(); ++i) {i->restore_prev_pos_range();}
	}
};

building_params_t global_building_params;

void buildings_file_err(string const &str, int &error) {
	cout << "Error reading buildings config option " << str << "." << endl;
	error = 1;
}
bool check_01(float v) {return (v >= 0.0 && v <= 1.0);}

bool parse_buildings_option(FILE *fp) {

	int error(0);
	char strc[MAX_CHARS] = {0};
	if (!read_str(fp, strc)) return 0;
	string const str(strc);

	// global parameters
	if (str == "flatten_mesh") {
		if (!read_bool(fp, global_building_params.flatten_mesh)) {buildings_file_err(str, error);}
	}
	else if (str == "num_place") {
		if (!read_uint(fp, global_building_params.num_place)) {buildings_file_err(str, error);}
	}
	else if (str == "num_tries") {
		if (!read_uint(fp, global_building_params.num_tries)) {buildings_file_err(str, error);}
	}
	else if (str == "ao_factor") {
		if (!read_zero_one_float(fp, global_building_params.ao_factor)) {buildings_file_err(str, error);}
	}
	else if (str == "sec_extra_spacing") {
		if (!read_float(fp, global_building_params.sec_extra_spacing)) {buildings_file_err(str, error);}
	}
	else if (str == "tt_only") {
		if (!read_bool(fp, global_building_params.tt_only)) {buildings_file_err(str, error);}
	}
	else if (str == "infinite_buildings") {
		if (!read_bool(fp, global_building_params.infinite_buildings)) {buildings_file_err(str, error);}
	}
	// material parameters
	else if (str == "range_translate") { // x,y only
		if (!(read_float(fp, global_building_params.range_translate.x) &&
			read_float(fp, global_building_params.range_translate.y))) {buildings_file_err(str, error);}
	}
	else if (str == "pos_range") {
		if (!read_cube(fp, global_building_params.cur_mat.pos_range, 1)) {buildings_file_err(str, error);}
	}
	else if (str == "place_radius") {
		if (!read_float(fp, global_building_params.cur_mat.place_radius)) {buildings_file_err(str, error);}
	}
	else if (str == "max_delta_z") {
		if (!read_float(fp, global_building_params.cur_mat.max_delta_z)) {buildings_file_err(str, error);}
	}
	else if (str == "min_level_height") {
		if (!read_float(fp, global_building_params.cur_mat.min_level_height)) {buildings_file_err(str, error);}
	}
	else if (str == "max_rot_angle") {
		if (!read_float(fp, global_building_params.cur_mat.max_rot_angle)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.max_rot_angle *= TO_RADIANS; // specified in degrees, stored in radians
	}
	else if (str == "split_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.split_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "cube_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.cube_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "round_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.round_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "alt_step_factor_prob") {
		if (!read_zero_one_float(fp, global_building_params.cur_mat.asf_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "min_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.min_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "max_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.max_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "min_sides") {
		if (!read_uint(fp, global_building_params.cur_mat.min_sides)) {buildings_file_err(str, error);}
		if (global_building_params.cur_mat.min_sides < 3) {buildings_file_err(str+" (< 3)", error);}
	}
	else if (str == "max_sides") {
		if (!read_uint(fp, global_building_params.cur_mat.max_sides)) {buildings_file_err(str, error);}
		if (global_building_params.cur_mat.max_sides < 3) {buildings_file_err(str+" (< 3)", error);}
	}
	else if (str == "min_flat_side_amt") {
		if (!read_float(fp, global_building_params.cur_mat.min_fsa)) {buildings_file_err(str, error);}
	}
	else if (str == "max_flat_side_amt") {
		if (!read_float(fp, global_building_params.cur_mat.max_fsa)) {buildings_file_err(str, error);}
	}
	else if (str == "min_alt_step_factor") {
		if (!read_float(fp, global_building_params.cur_mat.min_asf)) {buildings_file_err(str, error);}
	}
	else if (str == "max_alt_step_factor") {
		if (!read_float(fp, global_building_params.cur_mat.max_asf)) {buildings_file_err(str, error);}
	}
	else if (str == "size_range") {
		if (!read_cube(fp, global_building_params.cur_mat.sz_range)) {buildings_file_err(str, error);}
	}
	else if (str == "min_altitude") {
		if (!read_float(fp, global_building_params.cur_mat.min_alt)) {buildings_file_err(str, error);}
	}
	else if (str == "max_altitude") {
		if (!read_float(fp, global_building_params.cur_mat.max_alt)) {buildings_file_err(str, error);}
	}
	else if (str == "no_city") {
		if (!read_bool(fp, global_building_params.cur_mat.no_city)) {buildings_file_err(str, error);}
	}
	// material textures
	else if (str == "texture_mirror") {
		if (!read_bool(fp, global_building_params.tex_mirror)) {buildings_file_err(str, error);}
	}
	else if (str == "texture_inv_y") {
		if (!read_bool(fp, global_building_params.tex_inv_y)) {buildings_file_err(str, error);}
	}
	else if (str == "side_tscale") { // both X and Y
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale_x)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_tex.tscale_y = global_building_params.cur_mat.side_tex.tscale_x; // uniform
	}
	else if (str == "side_tscale_x") {
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale_x)) {buildings_file_err(str, error);}
	}
	else if (str == "side_tscale_y") {
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale_y)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_tscale") { // both X and Y
		if (!read_float(fp, global_building_params.cur_mat.roof_tex.tscale_x)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_tex.tscale_y = global_building_params.cur_mat.roof_tex.tscale_x; // uniform
	}
	else if (str == "side_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_tex.tid = get_texture_by_name(std::string(strc), 0, global_building_params.tex_inv_y, global_building_params.get_wrap_mir());
	}
	else if (str == "side_nm_tid") { // Warning: setting options such as tex_inv_y for textures that have already been loaded will have no effect!
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_tex.nm_tid = get_texture_by_name(std::string(strc), 1, global_building_params.tex_inv_y, global_building_params.get_wrap_mir());
	}
	else if (str == "roof_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_tex.tid = get_texture_by_name(std::string(strc), 0, global_building_params.tex_inv_y, global_building_params.get_wrap_mir());
	}
	else if (str == "roof_nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_tex.nm_tid = get_texture_by_name(std::string(strc), 1, global_building_params.tex_inv_y, global_building_params.get_wrap_mir());
	}
	// material colors
	else if (str == "side_color") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmin)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_color.cmax = global_building_params.cur_mat.side_color.cmin; // same
	}
	else if (str == "side_color_min") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmin)) {buildings_file_err(str, error);}
	}
	else if (str == "side_color_max") {
		if (!read_color(fp, global_building_params.cur_mat.side_color.cmax)) {buildings_file_err(str, error);}
	}
	else if (str == "side_color_grayscale_rand") {
		if (!read_float(fp, global_building_params.cur_mat.side_color.grayscale_rand)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmin)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_color.cmax = global_building_params.cur_mat.roof_color.cmin; // same
	}
	else if (str == "roof_color_min") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmin)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color_max") {
		if (!read_color(fp, global_building_params.cur_mat.roof_color.cmax)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_color_grayscale_rand") {
		if (!read_float(fp, global_building_params.cur_mat.roof_color.grayscale_rand)) {buildings_file_err(str, error);}
	}
	// windows
	else if (str == "window_width") {
		if (!read_float(fp, global_building_params.window_width) || !check_01(global_building_params.window_width)) {buildings_file_err(str, error);}
	}
	else if (str == "window_height") {
		if (!read_float(fp, global_building_params.window_height) || !check_01(global_building_params.window_height)) {buildings_file_err(str, error);}
	}
	else if (str == "window_xspace") {
		if (!read_float(fp, global_building_params.window_xspace) || !check_01(global_building_params.window_xspace)) {buildings_file_err(str, error);}
	}
	else if (str == "window_yspace") {
		if (!read_float(fp, global_building_params.window_yspace) || !check_01(global_building_params.window_yspace)) {buildings_file_err(str, error);}
	}
	else if (str == "window_xscale") {
		if (!read_float(fp, global_building_params.cur_mat.wind_xscale) || global_building_params.cur_mat.wind_xscale < 0.0) {buildings_file_err(str, error);}
	}
	else if (str == "window_yscale") {
		if (!read_float(fp, global_building_params.cur_mat.wind_yscale) || global_building_params.cur_mat.wind_yscale < 0.0) {buildings_file_err(str, error);}
	}
	else if (str == "window_xoff") {
		if (!read_float(fp, global_building_params.cur_mat.wind_xoff)) {buildings_file_err(str, error);}
	}
	else if (str == "window_yoff") {
		if (!read_float(fp, global_building_params.cur_mat.wind_yoff)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.wind_yoff *= -1.0; // invert Y
	}
	else if (str == "add_windows") { // per-material
		if (!read_bool(fp, global_building_params.cur_mat.add_windows)) {buildings_file_err(str, error);}
	}
	else if (str == "add_window_lights") { // per-material
		if (!read_bool(fp, global_building_params.cur_mat.add_wind_lights)) {buildings_file_err(str, error);}
	}
	else if (str == "house_prob") { // per-material
		if (!read_float(fp, global_building_params.cur_mat.house_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "house_scale_range") { // per-material
	if (!read_float(fp, global_building_params.cur_mat.house_scale_min) || !read_float(fp, global_building_params.cur_mat.house_scale_max)) {buildings_file_err(str, error);}
	}
	else if (str == "window_color") { // per-material
		if (!read_color(fp, global_building_params.cur_mat.window_color)) {buildings_file_err(str, error);}
	}
	// special commands
	else if (str == "probability") {
		if (!read_uint(fp, global_building_params.cur_prob)) {buildings_file_err(str, error);}
	}
	else if (str == "add_material") {global_building_params.add_cur_mat();}
	else {
		cout << "Unrecognized buildings keyword in input file: " << str << endl;
		error = 1;
	}
	return !error;
}


class building_draw_t;

struct building_geom_t { // describes the physical shape of a building
	unsigned num_sides;
	bool half_offset, is_pointed;
	bool door_dim, door_dir, door_part; // for door/window placement
	float rot_sin, rot_cos, flat_side_amt, alt_step_factor, start_angle; // rotation in XY plane, around Z (up) axis
	//float roof_recess;

	building_geom_t(unsigned ns=4, float rs=0.0, float rc=1.0) : num_sides(ns), half_offset(0), is_pointed(0), door_dim(0), door_dir(0), door_part(0),
		rot_sin(rs), rot_cos(rc), flat_side_amt(0.0), alt_step_factor(0.0), start_angle(0.0) {}
	bool is_rotated() const {return (rot_sin != 0.0);}
	bool is_cube()    const {return (num_sides == 4);}
	bool is_simple_cube()    const {return (is_cube() && !half_offset && flat_side_amt == 0.0 && alt_step_factor == 0.0);}
	bool use_cylinder_coll() const {return (num_sides > 8 && flat_side_amt == 0.0);} // use cylinder collision if not a cube, triangle, octagon, etc. (approximate)
};

struct tquad_with_ix_t : public tquad_t {
	enum {TYPE_ROOF=0, TYPE_WALL, TYPE_DOOR};
	unsigned type;
	tquad_with_ix_t(unsigned npts_=0) : tquad_t(npts_), type(0) {}
	tquad_with_ix_t(tquad_t const &t, unsigned type_) : tquad_t(t), type(type_) {}
};

struct building_t : public building_geom_t {

	unsigned mat_ix;
	bool is_house, has_antenna;
	colorRGBA side_color, roof_color, detail_color;
	cube_t bcube;
	vect_cube_t parts;
	vect_cube_t details; // cubes on the roof - antennas, AC units, etc.
	vector<tquad_with_ix_t> roof_tquads, doors;
	float ao_bcz2;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_), is_house(0), has_antenna(0), side_color(WHITE), roof_color(WHITE), detail_color(BLACK), ao_bcz2(0.0) {bcube.set_to_zeros();}
	bool is_valid() const {return !bcube.is_all_zeros();}
	colorRGBA get_avg_side_color  () const {return side_color.modulate_with(get_material().side_tex.get_avg_color());}
	colorRGBA get_avg_roof_color  () const {return roof_color.modulate_with(get_material().roof_tex.get_avg_color());}
	colorRGBA get_avg_detail_color() const {return detail_color.modulate_with(get_material().roof_tex.get_avg_color());}
	building_mat_t const &get_material() const {return global_building_params.get_material(mat_ix);}
	void gen_rotation(rand_gen_t &rgen);

	void set_z_range(float z1, float z2) {
		bcube.z1() = z1; bcube.z2() = z2;
		if (!parts.empty()) {parts[0].z1() = z1; parts[0].z2() = z2;}
	}
	bool check_part_contains_pt_xy(cube_t const &part, point const &pt, vector<point> &points) const;
	bool check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const {
		return (check_bcube_overlap_xy_one_dir(b, expand_rel, expand_abs, points) || b.check_bcube_overlap_xy_one_dir(*this, expand_rel, expand_abs, points));
	}
	bool check_sphere_coll(point const &pos, float radius, bool xy_only, vector<point> &points, vector3d *cnorm=nullptr) const {
		point pos2(pos);
		return check_sphere_coll(pos2, pos, zero_vector, radius, xy_only, points, cnorm);
	}
	bool check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector<point> &points, vector3d *cnorm=nullptr) const;
	unsigned check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t, vector<point> &points, bool occlusion_only=0, bool ret_any_pt=0) const;
	bool check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const;
	void calc_bcube_from_parts();
	void gen_geometry(int rseed1, int rseed2);
	void gen_house(cube_t const &base, rand_gen_t &rgen);
	void add_door(cube_t const &c, bool dim, bool dir);
	void gen_peaked_roof(cube_t const &top, float peak_height, bool dim);
	void gen_details(rand_gen_t &rgen);
	void gen_sloped_roof(rand_gen_t &rgen);
	void add_roof_to_bcube();
	void gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax);
	void get_all_drawn_verts(building_draw_t &bdraw) const;
	void get_all_drawn_window_verts(building_draw_t &bdraw, bool lights_pass) const;
private:
	bool check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const;
	void split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen);
	bool test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const;
};


void do_xy_rotate(float rot_sin, float rot_cos, point const &center, point &pos) {
	float const x(pos.x - center.x), y(pos.y - center.y); // translate to center
	pos.x = x*rot_cos - y*rot_sin + center.x;
	pos.y = y*rot_cos + x*rot_sin + center.y;
}
void do_xy_rotate_normal(float rot_sin, float rot_cos, point &n) {
	float const x(n.x), y(n.y);
	n.x = x*rot_cos - y*rot_sin;
	n.y = y*rot_cos + x*rot_sin;
}


class building_window_gen_t { // and doors?
	int window_tid, door_tid;
public:
	building_window_gen_t() : window_tid(-1), door_tid(-1) {}
	int get_window_tid() const {return window_tid;}
	
	int get_door_tid() {
		if (door_tid < 0) {door_tid = get_texture_by_name("textures/white_door.jpg");}
		if (door_tid < 0) {door_tid = WHITE_TEX;} // failed to load door texture - use a simple white texture
		return door_tid;
	}
	bool check_windows_texture() {
		if (!global_building_params.windows_enabled()) return 0;
		if (window_tid >= 0) return 1; // already generated
		gen_building_window_texture(global_building_params.get_window_width_fract(), global_building_params.get_window_height_fract());
		window_tid = BLDG_WINDOW_TEX;
		return 1;
	}
};
building_window_gen_t building_window_gen;


struct building_draw_utils {
	static void calc_normals(building_geom_t const &bg, vector<vector3d> &nv, unsigned ndiv) {
		assert(bg.flat_side_amt >= 0.0 && bg.flat_side_amt < 0.5); // generates a flat side
		assert(bg.alt_step_factor >= 0.0 && bg.alt_step_factor < 1.0);
		if (bg.flat_side_amt > 0.0) {assert(ndiv > 4);} // should be at least 5 sides, 6-8 is better
		float const ndiv_inv(1.0/ndiv), css(TWO_PI*ndiv_inv*(1.0f - bg.flat_side_amt));
		float sin_ds[2], cos_ds[2];

		if (bg.alt_step_factor > 0.0) { // alternate between large and small steps (cube with angled corners, etc.)
			assert(!(ndiv&1));
			float const css_v[2] = {css*(1.0f + bg.alt_step_factor), css*(1.0f - bg.alt_step_factor)};
			UNROLL_2X(sin_ds[i_] = sin(css_v[i_]); cos_ds[i_] = cos(css_v[i_]);)
		}
		else { // uniform side length
			sin_ds[0] = sin_ds[1] = sin(css);
			cos_ds[0] = cos_ds[1] = cos(css);
		}
		float sin_s(0.0), cos_s(1.0), angle0(bg.start_angle); // start at 0.0
		if (bg.half_offset) {angle0 = 0.5*css;} // for cube
		if (angle0 != 0.0) {sin_s = sin(angle0); cos_s = cos(angle0);} // uncommon case
		nv.resize(ndiv);

		for (unsigned S = 0; S < ndiv; ++S) { // build normals table
			bool const d(S&1);
			float const s(sin_s), c(cos_s);
			nv[S].assign(s, c, 0.0);
			sin_s = s*cos_ds[d] + c*sin_ds[d];
			cos_s = c*cos_ds[d] - s*sin_ds[d];
		}
	}

	static void calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, vector<point> &pts, float expand=0.0) {
		calc_normals(bg, pts, bg.num_sides);
		vector3d const sz(bcube.get_size());
		point const cc(bcube.get_cube_center());
		float const rx(0.5*sz.x + expand), ry(0.5*sz.y + expand); // expand polygon by sphere radius
		for (unsigned i = 0; i < bg.num_sides; ++i) {pts[i].assign((cc.x + rx*pts[i].x), (cc.y + ry*pts[i].y), 0.0);} // convert normals to points
	}
};


#define EMIT_VERTEX() \
	vert.v = pt*sz + llc; \
	vert.t[ st] = tscale[ st]*(vert.v[d] + tex_vert_off[d]); \
	vert.t[!st] = tscale[!st]*(vert.v[i] + tex_vert_off[i]); \
	vert.t[0] += tex.txoff; \
	vert.t[1] += tex.tyoff; \
	if (apply_ao) {vert.copy_color(cw[pt.z == 1]);} \
	if (bg.is_rotated()) {do_xy_rotate(bg.rot_sin, bg.rot_cos, center, vert.v);} \
	verts.push_back(vert);

typedef vector<vert_norm_comp_tc_color> vect_vnctcc_t;

class building_draw_t {

	class draw_block_t {
		struct vert_ix_pair {
			unsigned qix, tix; // {quads, tris}
			vert_ix_pair(unsigned qix_, unsigned tix_) : qix(qix_), tix(tix_) {}
			bool operator==(vert_ix_pair const &v) const {return (qix == v.qix && tix == v.tix);}
		};
		// Note: not easy to use vao_manager_t due to upload done before active shader + shadow vs. geometry pass, but we can use vao_wrap_t's directly
		vbo_wrap_t qvbo, tvbo;
		vao_wrap_t qvao, tvao, sqvao, stvao; // regular + shadow
		vector<vert_ix_pair> pos_by_tile; // {quads, tris}
	public:
		tid_nm_pair_t tex;
		vect_vnctcc_t quad_verts, tri_verts;

		void draw_geom_range(bool shadow_only, vert_ix_pair const &vstart, vert_ix_pair const &vend) { // use VBO rendering
			if (vstart == vend) return; // empty range - no verts for this tile
			if (!shadow_only) {tex.set_gl();}
			ensure_vbos();

			if (qvbo.vbo_valid() && vstart.qix != vend.qix) {
				(shadow_only ? sqvao : qvao).create_from_vbo<vert_norm_comp_tc_color>(qvbo, 1, 1); // setup_pointers=1, always_bind=1
				assert(vstart.qix < vend.qix);
				draw_quads_as_tris((vend.qix - vstart.qix), vstart.qix);
			}
			if (tvbo.vbo_valid() && vstart.tix != vend.tix) {
				(shadow_only ? stvao : tvao).create_from_vbo<vert_norm_comp_tc_color>(tvbo, 1, 1); // setup_pointers=1, always_bind=1
				assert(vstart.tix < vend.tix);
				glDrawArrays(GL_TRIANGLES, vstart.tix, (vend.tix - vstart.tix));
			}
			vao_manager_t::post_render();
		}
		void draw_all_geom(bool shadow_only) {
			if (pos_by_tile.empty()) return; // nothing to draw for this block/texture
			draw_geom_range(shadow_only, pos_by_tile.front(), pos_by_tile.back());
		}
		void draw_geom_tile(unsigned tile_id) {
			if (pos_by_tile.empty()) return; // nothing to draw for this block/texture
			assert(tile_id+1 < pos_by_tile.size()); // tile and next tile must be valid indices
			draw_geom_range(0, pos_by_tile[tile_id], pos_by_tile[tile_id+1]); // shadow_only=0
		}
		void ensure_vbos() {
			assert((quad_verts.size()%4) == 0);
			assert((tri_verts.size()%3) == 0);
			if (!quad_verts.empty()) {qvbo.create_and_upload(quad_verts, 0, 1);}
			if (! tri_verts.empty()) {tvbo.create_and_upload( tri_verts, 0, 1);}
		}
		void upload_to_vbos() {
			ensure_vbos();

			if (no_store_model_textures_in_memory) { // no longer needed, unless VBO needs to be recreated for maximize/F9; use same condition as model textures
				clear_cont(quad_verts);
				clear_cont(tri_verts);
			}
		}
		void register_tile_id(unsigned tid) {
			if (tid+1 == pos_by_tile.size()) return; // already saw this tile
			assert(tid >= pos_by_tile.size()); // tid must be strictly increasing
			pos_by_tile.resize(tid+1, vert_ix_pair(quad_verts.size(), tri_verts.size())); // push start of new range back onto all previous tile slots
		}
		void finalize(unsigned num_tiles) {
			if (empty()) return; // nothing to do
			register_tile_id(num_tiles); // add terminator
			remove_excess_cap(pos_by_tile);

			if (!no_store_model_textures_in_memory) {
				remove_excess_cap(quad_verts);
				remove_excess_cap(tri_verts);
			}
		}
		void clear_verts() {quad_verts.clear(); tri_verts.clear(); pos_by_tile.clear();}
		void clear_vbos () {qvbo.clear(); tvbo.clear(); qvao.clear(); tvao.clear(); sqvao.clear(); stvao.clear();}
		void clear() {clear_vbos(); clear_verts();}
		bool empty() const {return (quad_verts.empty() && tri_verts.empty());}
		unsigned num_verts() const {return (quad_verts.size() + tri_verts.size());}
		unsigned num_tris () const {return (quad_verts.size()/2 + tri_verts.size()/3);} // Note: 1 quad = 4 verts = 2 triangles
	}; // end draw_block_t
	vector<draw_block_t> to_draw; // one per texture, assumes tids are dense

	vect_vnctcc_t &get_verts(tid_nm_pair_t const &tex, bool quads_or_tris=0) { // default is quads
		unsigned const ix((tex.tid >= 0) ? (tex.tid+1) : 0);
		if (ix >= to_draw.size()) {to_draw.resize(ix+1);}
		draw_block_t &block(to_draw[ix]);
		block.register_tile_id(cur_tile_id);
		if (block.empty()) {block.tex = tex;} // copy material first time
		else {assert(block.tex.nm_tid == tex.nm_tid);} // else normal maps must agree
		return (quads_or_tris ? block.tri_verts : block.quad_verts);
	}
	static void setup_ao_color(colorRGBA const &color, float bcz1, float ao_bcz2, float z1, float z2, color_wrapper cw[2], vert_norm_comp_tc_color &vert, bool no_ao) {
		if (!no_ao && global_building_params.ao_factor > 0.0) {
			float const dz_mult(global_building_params.ao_factor/(ao_bcz2 - bcz1));
			UNROLL_2X(cw[i_].set_c4(color*((1.0 - global_building_params.ao_factor) + dz_mult*((i_ ? z2 : z1) - bcz1)));)
		} else {vert.set_c4(color);} // color is shared across all verts
	}
	vector<vector3d> normals; // reused across add_cylinder() calls
	point cur_camera_pos;

public:
	unsigned cur_tile_id;
	building_draw_t() : cur_camera_pos(zero_vector), cur_tile_id(0) {}
	void init_draw_frame() {cur_camera_pos = get_camera_pos();} // capture camera pos during non-shadow pass to use for shadow pass
	bool empty() const {return to_draw.empty();}

	void add_cylinder(building_geom_t const &bg, point const &pos, point const &rot_center, float height, float rx, float ry, float bcz1, float ao_bcz2,
		tid_nm_pair_t const &tex, colorRGBA const &color, unsigned dim_mask, bool no_ao, bool clip_windows)
	{
		unsigned const ndiv(bg.num_sides); // Note: no LOD
		assert(ndiv >= 3);
		bool const smooth_normals(ndiv >= 16); // cylinder vs. N-gon
		float const z_top(pos.z + height), tscale_x(2.0*tex.tscale_x), tscale_y(2.0*tex.tscale_y); // adjust for local vs. global space change
		bool const apply_ao(!no_ao && global_building_params.ao_factor > 0.0);
		vert_norm_comp_tc_color vert;
		color_wrapper cw[2];
		setup_ao_color(color, bcz1, ao_bcz2, pos.z, z_top, cw, vert, no_ao);
		float tex_pos[2] = {0.0, 1.0};
		building_draw_utils::calc_normals(bg, normals, ndiv);
		UNROLL_2X(tex_pos[i_] = ((i_ ? z_top : pos.z) - bcz1);)

		if (dim_mask & 3) { // draw sides
			auto &verts(get_verts(tex)); // Note: cubes are drawn with quads, so we want to emit quads here
			float tot_perim(0.0), cur_perim[2] = {0.0, 0.0};
			for (unsigned S = 0; S < ndiv; ++S) {tot_perim += p2p_dist(normals[S], normals[(S+1)%ndiv]);}
			float const tscale_mult(TWO_PI*sqrt((rx*rx + ry*ry)/2.0)/tot_perim);
				
			for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data quads
				vector3d const &n1(normals[S]), &n2(normals[(S+1)%ndiv]);
				cur_perim[0]  = cur_perim[1];
				cur_perim[1] += p2p_dist(n1, n2);
				vector3d normal(n1 + n2); normal.x *= ry; normal.y *= rx; // average the two vertex normals for the flat face normal
				if (bg.is_rotated()) {do_xy_rotate_normal(bg.rot_sin, bg.rot_cos, normal);}
				bool const cur_smooth_normals(smooth_normals && (bg.flat_side_amt == 0.0 || S+1 != ndiv)); // flat side of cylindrical building is not smooth
				if (!cur_smooth_normals) {vert.set_norm(normal.get_norm());}

				for (unsigned d = 0; d < 2; ++d) {
					vector3d const &n(d ? n2 : n1);
					vert.t[0] = tscale_x*cur_perim[d]*tscale_mult + tex.txoff; // Note: could try harder to ensure an integer multiple to fix seams, but not a problem in practice
					
					if (cur_smooth_normals) {
						vector3d normal(n); normal.x *= ry; normal.y *= rx; // scale normal by radius (swapped)
						if (bg.is_rotated()) {do_xy_rotate_normal(bg.rot_sin, bg.rot_cos, normal);}
						vert.set_norm(normal.get_norm());
					}
					vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), 0.0);
					if (bg.is_rotated()) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, vert.v);}

					for (unsigned e = 0; e < 2; ++e) {
						vert.v.z = ((d^e) ? z_top : pos.z);
						vert.t[1] = tscale_y*tex_pos[d^e] + tex.tyoff;
						if (apply_ao) {vert.copy_color(cw[d^e]);}
						verts.push_back(vert);
					}
					if (clip_windows) {clip_low_high(verts[verts.size()-1].t[1], verts[verts.size()-2].t[1]);} // is this necessary?
				} // for d
			} // for S
		} // end draw sides
		if (dim_mask & 4) { // draw end(s) / roof
			auto &tri_verts(get_verts(tex, 1));
			
			for (unsigned d = 0; d < 2; ++d) { // bottom, top
				vert.set_norm(d ? plus_z : -plus_z);
				if (apply_ao) {vert.copy_color(cw[d]);}
				vert_norm_comp_tc_color center(vert);
				center.t[0] = center.t[1] = 0.0; // center of texture space for this disk
				center.v = pos;
				if (d) {center.v.z += height;}
				if (bg.is_rotated()) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, center.v);}

				for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data triangles
					tri_verts.push_back(center);

					for (unsigned e = 0; e < 2; ++e) {
						if (S > 0 && e == 0) {tri_verts.push_back(tri_verts[tri_verts.size()-2]); continue;} // reuse prev vertex
						vector3d const &n(normals[(S+e)%ndiv]);
						vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), center.v.z);
						vert.t[0] = tscale_x*n[0]; vert.t[1] = tscale_y*n[1];
						if (bg.is_rotated()) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, vert.v);}
						tri_verts.push_back(vert);
					}
				} // for S
			} // for d
		} // end draw end(s)
	}

	void add_tquad(building_geom_t const &bg, tquad_with_ix_t const &tquad, cube_t const &bcube, tid_nm_pair_t const &tex, colorRGBA const &color) {
		assert(tquad.npts == 3 || tquad.npts == 4); // triangles or quads
		auto &verts(get_verts(tex, (tquad.npts == 3))); // 0=quads, 1=tris
		point const center(!bg.is_rotated() ? all_zeros : bcube.get_cube_center()); // rotate about bounding cube / building center
		vert_norm_comp_tc_color vert;
		float tsx(0.0), tsy(0.0), tex_off(0.0);
		bool dim(0);

		if (tquad.type == tquad_with_ix_t::TYPE_WALL) { // side/wall
			tsx = 2.0f*tex.tscale_x; tsy = 2.0f*tex.tscale_y; // adjust for local vs. global space change
			dim = (tquad.pts[0].x == tquad.pts[1].x);
			if (world_mode != WMODE_INF_TERRAIN) {tex_off = (dim ? yoff2*DY_VAL : xoff2*DX_VAL);}
		}
		else if (tquad.type == tquad_with_ix_t::TYPE_ROOF) { // roof
			float const denom(0.5*(bcube.get_dx() + bcube.get_dy()));
			tsx = tex.tscale_x/denom; tsy = tex.tscale_y/denom;
		}
		vert.set_c4(color);
		vector3d normal(tquad.get_norm());
		if (bg.is_rotated()) {do_xy_rotate_normal(bg.rot_sin, bg.rot_cos, normal);}
		vert.set_norm(normal);
		
		for (unsigned i = 0; i < tquad.npts; ++i) {
			vert.v = tquad.pts[i];

			if (tquad.type == tquad_with_ix_t::TYPE_WALL) { // side/wall
				vert.t[0] = (vert.v[dim] + tex_off)*tsx; // use nonzero width dim
				vert.t[1] = vert.v.z*tsy;
			}
			else if (tquad.type == tquad_with_ix_t::TYPE_ROOF) { // roof
				vert.t[0] = (vert.v.x - bcube.x1())*tsx; // varies from 0.0 and bcube x1 to 1.0 and bcube x2
				vert.t[1] = (vert.v.y - bcube.y1())*tsy; // varies from 0.0 and bcube y1 to 1.0 and bcube y2
			}
			else if (tquad.type == tquad_with_ix_t::TYPE_DOOR) { // door - textured from (0,0) to (1,1)
				vert.t[0] = float(i == 1 || i == 2);
				vert.t[1] = float(i == 2 || i == 3);
			}
			else {assert(0);}
			if (bg.is_rotated()) {do_xy_rotate(bg.rot_sin, bg.rot_cos, center, vert.v);}
			verts.push_back(vert);
		}
	}

	void add_section(building_geom_t const &bg, cube_t const &cube, cube_t const &bcube, float ao_bcz2, tid_nm_pair_t const &tex,
		colorRGBA const &color, unsigned dim_mask, bool skip_bottom, bool skip_top, bool no_ao, bool clip_windows, float door_ztop=0.0)
	{
		assert(bg.num_sides >= 3); // must be nonzero volume
		point const center(!bg.is_rotated() ? all_zeros : bcube.get_cube_center()); // rotate about bounding cube / building center
		vector3d const sz(cube.get_size());

		if (bg.num_sides != 4) { // not a cube, use cylinder
			assert(door_ztop == 0.0); // not supported
			point const ccenter(cube.get_cube_center()), pos(ccenter.x, ccenter.y, cube.z1());
			//float const rscale(0.5*((num_sides <= 8) ? SQRT2 : 1.0)); // larger for triangles/cubes/hexagons/octagons (to ensure overlap/connectivity), smaller for cylinders
			float const rscale(0.5); // use shape contained in bcube so that bcube tests are correct, since we're not creating L/T/U shapes for this case
			add_cylinder(bg, pos, center, sz.z, rscale*sz.x, rscale*sz.y, bcube.z1(), ao_bcz2, tex, color, dim_mask, no_ao, clip_windows);
			return;
		}
		// else draw as a cube (optimized flow)
		auto &verts(get_verts(tex, bg.is_pointed)); // bg.is_pointed ? tris : quads
		vector3d const llc(cube.get_llc()); // move origin from center to min corner
		vert_norm_comp_tc_color vert;
		if (bg.is_pointed) {dim_mask &= 3;} // mask off z-dim since pointed objects (antenna) have no horizontal surfaces
		float const tscale[2] = {2.0f*tex.tscale_x, 2.0f*tex.tscale_y}; // adjust for local vs. global space change
		bool const apply_ao(!no_ao && global_building_params.ao_factor > 0.0);
		color_wrapper cw[2];
		setup_ao_color(color, bcube.z1(), ao_bcz2, cube.d[2][0], cube.d[2][1], cw, vert, no_ao);
		vector3d const tex_vert_off(((world_mode == WMODE_INF_TERRAIN) ? zero_vector : vector3d(xoff2*DX_VAL, yoff2*DY_VAL, 0.0)));
		
		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const n((i+2)%3), d((i+1)%3), st(i&1); // direction of normal
			if (!(dim_mask & (1<<n))) continue;

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (skip_bottom && n == 2 && j == 0) continue; // skip bottom side
				if (skip_top    && n == 2 && j == 1) continue; // skip top    side
				
				if (n < 2 && bg.is_rotated()) { // XY only
					vector3d norm; norm.z = 0.0;
					if (n == 0) {norm.x =  bg.rot_cos; norm.y = bg.rot_sin;} // X
					else        {norm.x = -bg.rot_sin; norm.y = bg.rot_cos;} // Y
					vert.set_norm(j ? norm : -norm);
				}
				else {
					vert.n[i] = 0; vert.n[d] = 0; vert.n[n] = (j ? 127 : -128); // -1.0 or 1.0
				}
				unsigned const ix(verts.size()); // first vertex of this quad/triangle
				point pt; pt[n] = j; // in direction of normal

				if (bg.is_pointed) {
					assert(door_ztop == 0.0); // not supported
					pt[!n] = !j; pt.z = 0;
					EMIT_VERTEX(); // bottom low
					pt[!n] = j;
					EMIT_VERTEX(); // bottom high
					pt[!n] = 0.5; pt[n] = 0.5; pt.z = 1;
					EMIT_VERTEX(); // top
					vector3d normal;
					get_normal(verts[ix+0].v, verts[ix+1].v, verts[ix+2].v, normal, 1); // update with correct normal
					vert.set_norm(n ? -normal : normal);
					UNROLL_3X(verts[ix+i_].set_norm(vert);)
					continue; // no windows/clipping
				}
				pt[d] = 0;
				pt[i] = !j; // need to orient the vertices differently for each side
				//if (bg.roof_recess > 0.0 && n == 2 && j == 1) {pt.z -= bg.roof_recess*cube.get_dz();}
				EMIT_VERTEX(); // 0 !j
				pt[i] = j;
				EMIT_VERTEX(); // 0 j
				pt[d] = 1;
				EMIT_VERTEX(); // 1 j
				pt[i] = !j;
				EMIT_VERTEX(); // 1 !j

				if (door_ztop != 0.0 && n == unsigned(bg.door_dim) && j == unsigned(bg.door_dir)) { // not the cleanest way to do this ...
					float const door_height(door_ztop - cube.z1()), offset(0.03*(bg.door_dir ? 1.0 : -1.0)*door_height);

					for (unsigned k = ix; k < ix+4; ++k) {
						auto &v(verts[k]);
						float const delta(door_ztop - v.v.z);
						max_eq(v.v.z, door_ztop); // make all windows start above the door
						v.v[bg.door_dim] += offset; // move slightly away from the house wall to avoid z-fighting (vertex is different from building and won't have same depth)
						if (delta > 0.0) {v.t[1] += tscale[1]*delta;} // recalculate tex coord
					}
				}
				if (clip_windows && n < 2) { // clip the quad that was just added (side of building)
					clip_low_high(verts[ix+0].t[!st], verts[ix+1].t[!st]);
					clip_low_high(verts[ix+2].t[!st], verts[ix+3].t[!st]);
					clip_low_high(verts[ix+0].t[ st], verts[ix+3].t[ st]);
					clip_low_high(verts[ix+1].t[ st], verts[ix+2].t[ st]);
				}
			} // for j
		} // for i
	}

	static void clip_low_high(float &t0, float &t1) {
		if (fabs(t0 - t1) < 0.5) {t0 = t1 = 0.0;} // too small to have a window
		else {t0 = round_fp(t0); t1 = round_fp(t1);} // Note: round() is much faster than nearbyint(), and round_fp() is faster than round()
	}

	unsigned num_verts() const {
		unsigned num(0);
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {num += i->num_verts();}
		return num;
	}
	unsigned num_tris() const {
		unsigned num(0);
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {num += i->num_tris();}
		return num;
	}
	void upload_to_vbos() {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->upload_to_vbos();}}
	void clear_vbos    () {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->clear_vbos();}}
	void clear         () {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->clear();}}
	unsigned get_num_draw_blocks() const {return to_draw.size();}
	void finalize(unsigned num_tiles) {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->finalize(num_tiles);}}
	void draw          (bool shadow_only) {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->draw_all_geom (shadow_only);}}
	void draw_tile     (unsigned tile_id) {for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->draw_geom_tile(tile_id);}}
	void draw_block(unsigned ix, bool shadow_only) {if (ix < to_draw.size()) {to_draw[ix].draw_all_geom(shadow_only);}}
};


void building_t::split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen) {

	// generate L, T, U, or H shape
	point const llc(seed_cube.get_llc()), sz(seed_cube.get_size());
	int const shape(rand()%8); // 0-7
	bool const is_h(shape >= 7);
	bool const dim(rgen.rand_bool()); // {x,y}
	bool const dir(is_h ? 1 : rgen.rand_bool()); // {neg,pos} - H-shape is always pos
	float const div(is_h ? rgen.rand_uniform(0.2, 0.4) : rgen.rand_uniform(0.3, 0.7)), s1(rgen.rand_uniform(0.2, 0.4)), s2(rgen.rand_uniform(0.6, 0.8)); // split pos in 0-1 range
	float const dpos(llc[dim] + div*sz[dim]), spos1(llc[!dim] + s1*sz[!dim]), spos2(llc[!dim] + s2*sz[!dim]); // split pos in cube space
	unsigned const start(parts.size()), num((shape >= 6) ? 3 : 2);
	parts.resize(start+num, seed_cube);
	parts[start+0].d[dim][ dir] = dpos; // full width part
	parts[start+1].d[dim][!dir] = dpos; // partial width part

	switch (shape) {
	case 0: case 1: case 2: case 3: // L
		parts[start+1].d[!dim][shape>>1] = ((shape&1) ? spos2 : spos1);
		break;
	case 4: case 5: // T
		parts[start+1].d[!dim][0] = spos1;
		parts[start+1].d[!dim][1] = spos2;
		break;
	case 6: // U
		parts[start+2].d[ dim][!dir] = dpos; // partial width part
		parts[start+1].d[!dim][1  ] = spos1;
		parts[start+2].d[!dim][0  ] = spos2;
		break;
	case 7: { // H
		float const dpos2(llc[dim] + (1.0 - div)*sz[dim]); // other end
		parts[start+1].d[ dim][ dir] = dpos2;
		parts[start+1].d[!dim][ 0  ] = spos1;
		parts[start+1].d[!dim][ 1  ] = spos2;
		parts[start+2].d[ dim][!dir] = dpos2; // full width part
		break;
	}
	default: assert(0);
	}
}

void building_t::gen_rotation(rand_gen_t &rgen) {

	float const max_rot_angle(get_material().max_rot_angle);
	if (max_rot_angle == 0.0) return;
	float const rot_angle(rgen.rand_uniform(0.0, max_rot_angle));
	rot_sin = sin(rot_angle);
	rot_cos = cos(rot_angle);
	parts.clear();
	parts.push_back(bcube); // this is the actual building base
	cube_t const &bc(parts.back());
	point const center(bc.get_cube_center());

	for (unsigned i = 0; i < 4; ++i) {
		point corner(bc.d[0][i&1], bc.d[1][i>>1], bc.d[2][i&1]);
		do_xy_rotate(rot_sin, rot_cos, center, corner);
		if (i == 0) {bcube.set_from_point(corner);} else {bcube.union_with_pt(corner);} // Note: detail cubes are excluded
	}
}

// Note: only checks for point (x,y) value contained in one cube/N-gon/cylinder; assumes pt has already been rotated into local coordinate frame
bool building_t::check_part_contains_pt_xy(cube_t const &part, point const &pt, vector<point> &points) const {

	if (!part.contains_pt_xy(pt)) return 0; // check bounding cube
	if (is_simple_cube()) return 1; // that's it
	building_draw_utils::calc_poly_pts(*this, part, points);
	return point_in_polygon_2d(pt.x, pt.y, &points.front(), points.size(), 0, 1); // 2D x/y containment
}

bool building_t::check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const { // can be called before levels/splits are created

	if (expand_rel == 0.0 && expand_abs == 0.0 && !bcube.intersects(b.bcube)) return 0;
	if (!is_rotated() && !b.is_rotated()) return 1; // above check is exact, top-level bcube check up to the caller
	point const center1(b.bcube.get_cube_center()), center2(bcube.get_cube_center());
	//if (b.bcube.contains_pt_xy(center2) || bcube.contains_pt_xy(center1)) return 1; // slightly faster?
	
	for (auto p1 = b.parts.begin(); p1 != b.parts.end(); ++p1) {
		point pts[9]; // {center, 00, 10, 01, 11, x0, x1, y0, y1}
		cube_t c_exp, c_exp_rot;

		if (b.parts.size() == 1) { // single cube
			// don't need to calculate rotation because we know we're rotating about its center and bcube is correct; results are slightly different from the code below
			pts[0] = center1; // no need to rotate
			c_exp = b.bcube; // no need to rotate
			c_exp.expand_by_xy(expand_rel*p1->get_size() + vector3d(expand_abs, expand_abs, expand_abs));
			for (unsigned i = 0; i < 4; ++i) {pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0);} // XY only
		}
		else {
			pts[0] = p1->get_cube_center();
			do_xy_rotate(b.rot_sin, b.rot_cos, center1, pts[0]); // rotate into global space
			c_exp = cube_t(*p1);
			c_exp.expand_by_xy(expand_rel*p1->get_size() + vector3d(expand_abs, expand_abs, expand_abs));

			for (unsigned i = 0; i < 4; ++i) { // {00, 10, 01, 11}
				pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0); // XY only
				do_xy_rotate(b.rot_sin, b.rot_cos, center1, pts[i+1]); // rotate into global space
			}
		}
		for (unsigned i = 0; i < 5; ++i) {do_xy_rotate(-rot_sin, rot_cos, center2, pts[i]);} // inverse rotate into local coord space - negate the sine term
		c_exp_rot.set_from_points(pts+1, 4); // use points 1-4
		pts[5] = 0.5*(pts[1] + pts[3]); // x0 edge center
		pts[6] = 0.5*(pts[2] + pts[4]); // x1 edge center
		pts[7] = 0.5*(pts[1] + pts[2]); // y0 edge center
		pts[8] = 0.5*(pts[3] + pts[4]); // y1 edge center
		
		for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) {
			if (c_exp_rot.contains_pt_xy(p2->get_cube_center())) return 1; // quick and easy test for heavy overlap

			for (unsigned i = 0; i < 9; ++i) {
				if (check_part_contains_pt_xy(*p2, pts[i], points)) return 1; // Note: building geometry is likely not yet generated, this check should be sufficient
				//if (p2->contains_pt_xy(pts[i])) return 1;
			}
		}
	} // for p1
	return 0;
}

bool building_t::test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const {

	building_draw_utils::calc_poly_pts(*this, part, points); // without the expand
	point quad_pts[4]; // quads
	bool updated(0);

	// FIXME: if the player is moving too quickly, the intersection with a side polygon may be missed,
	// which allows the player to travel through the building, but using a line intersection test from p_past2 to pos has other problems
	for (unsigned S = 0; S < num_sides; ++S) { // generate vertex data quads
		for (unsigned d = 0, ix = 0; d < 2; ++d) {
			point const &p(points[(S+d)%num_sides]);
			for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, part.d[2][d^e]);}
		}
		vector3d const normal(get_poly_norm(quad_pts));
		float const rdist(dot_product_ptv(normal, pos, quad_pts[0]));
		if (rdist < 0.0 || rdist >= radius) continue; // too far or wrong side
		if (!sphere_poly_intersect(quad_pts, 4, pos, normal, rdist, radius)) continue;
		pos += normal*(radius - rdist);
		if (cnorm) {*cnorm = normal;}
		updated = 1;
	} // for S
	if (updated) return 1;
	
	if (max(pos.z, p_last.z) > part.z2() && point_in_polygon_2d(pos.x, pos.y, &points.front(), num_sides, 0, 1)) { // test top plane (sphere on top of polygon?)
		pos.z = part.z2() + radius; // make sure it doesn't intersect the roof
		if (cnorm) {*cnorm = plus_z;}
		return 1;
	}
	return 0;
}

bool building_t::check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only, vector<point> &points, vector3d *cnorm_ptr) const {

	if (!is_valid()) return 0; // invalid building
	point p_int;
	vector3d cnorm; // unused
	unsigned cdir(0); // unused
	if (!sphere_cube_intersect(pos, radius, (bcube + xlate), p_last, p_int, cnorm, cdir, 1, xy_only)) return 0;
	point pos2(pos), p_last2(p_last), center;
	bool had_coll(0);
	
	if (is_rotated()) {
		center = bcube.get_cube_center() + xlate;
		do_xy_rotate(-rot_sin, rot_cos, center, pos2); // inverse rotate - negate the sine term
		do_xy_rotate(-rot_sin, rot_cos, center, p_last2);
	}
	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (xy_only && i->d[2][0] > bcube.d[2][0]) break; // only need to check first level in this mode
		if (!xy_only && ((pos2.z + radius < i->d[2][0] + xlate.z) || (pos2.z - radius > i->d[2][1] + xlate.z))) continue; // test z overlap

		if (use_cylinder_coll()) {
			point const cc(i->get_cube_center() + xlate);
			float const crx(0.5*i->get_dx()), cry(0.5*i->get_dy()), r_sum(radius + max(crx, cry));
			if (!dist_xy_less_than(pos2, cc, r_sum)) continue; // no intersection

			if (fabs(crx - cry) < radius) { // close to a circle
				if (p_last2.z > i->d[2][1] + xlate.z && dist_xy_less_than(pos2, cc, max(crx, cry))) {
					pos2.z = i->z2() + radius; // make sure it doesn't intersect the roof
					if (cnorm_ptr) {*cnorm_ptr = plus_z;}
				}
				else { // side coll
					vector2d const d((pos2.x - cc.x), (pos2.y - cc.y));
					float const mult(r_sum/d.mag());
					pos2.x = cc.x + mult*d.x;
					pos2.y = cc.y + mult*d.y;
					if (cnorm_ptr) {*cnorm_ptr = vector3d(d.x, d.y, 0.0).get_norm();} // no z-component
				}
				had_coll = 1;
			}
			else {
				had_coll |= test_coll_with_sides(pos2, p_last2, radius, (*i + xlate), points, cnorm_ptr); // use polygon collision test
			}
		}
		else if (num_sides != 4) { // triangle, hexagon, octagon, etc.
			had_coll |= test_coll_with_sides(pos2, p_last2, radius, (*i + xlate), points, cnorm_ptr);
		}
		else if (sphere_cube_int_update_pos(pos2, radius, (*i + xlate), p_last2, 1, xy_only, cnorm_ptr)) { // cube
			had_coll = 1; // flag as colliding, continue to look for more collisions (inside corners)
		}
	} // for i
	for (auto i = details.begin(); i != details.end(); ++i) {
		if (sphere_cube_int_update_pos(pos2, radius, (*i + xlate), p_last2, 1, xy_only, cnorm_ptr)) {had_coll = 1;} // cube, flag as colliding
	}
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) { // Note: doesn't really work with a pointed roof
		point const pos_xlate(pos2 - xlate);
		vector3d const normal(i->get_norm());
		float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

		if (fabs(rdist) < radius && sphere_poly_intersect(i->pts, i->npts, pos_xlate, normal, rdist, radius)) {
			pos2 += normal*(radius - rdist); // update current pos
			had_coll = 1; // flag as colliding
			if (cnorm_ptr) {*cnorm_ptr = ((normal.z < 0.0) ? -1.0 : 1.0)*normal;} // make sure normal points up
			break; // only use first colliding tquad
		}
	}
	if (!had_coll) return 0; // Note: no collisions with windows or doors, since they're colinear with walls

	if (is_rotated()) {
		do_xy_rotate(rot_sin, rot_cos, center, pos2); // rotate back around center
		if (cnorm_ptr) {do_xy_rotate(rot_sin, rot_cos, all_zeros, *cnorm_ptr);} // rotate back (pure rotation)
	}
	pos = pos2;
	return had_coll;
}

unsigned building_t::check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t, vector<point> &points, bool occlusion_only, bool ret_any_pt) const {

	if (!check_line_clip(p1-xlate, p2-xlate, bcube.d)) return 0; // no intersection
	point p1r(p1), p2r(p2);
	float tmin(0.0), tmax(1.0);
	unsigned coll(0); // 0=none, 1=side, 2=roof, 3=details

	if (is_rotated()) {
		point const center(bcube.get_cube_center() + xlate);
		do_xy_rotate(-rot_sin, rot_cos, center, p1r); // inverse rotate - negate the sine term
		do_xy_rotate(-rot_sin, rot_cos, center, p2r);
	}
	p1r -= xlate; p2r -= xlate;
	float const pzmin(min(p1r.z, p2r.z)), pzmax(max(p1r.z, p2r.z));
	bool const vert(p1r.x == p2r.x && p1r.y == p2r.y);

	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (pzmin > i->z2() || pzmax < i->z1()) continue; // no overlap in z
		bool hit(0);

		if (use_cylinder_coll()) { // vertical cylinder
			// Note: we know the line intersects the cylinder's bcube, and there's a good chance it intersects the cylinder, so we don't need any expensive early termination cases here
			point const cc(i->get_cube_center());
			vector3d const csz(i->get_size());
			
			if (vert) { // vertical line + vertical cylinder optimization + handling of ellipsoids
				if (!point_in_ellipse(p1r, cc, 0.5*csz.x, 0.5*csz.y)) continue; // no intersection (below test should return true as well)
				tmin = (i->z2() - p1r.z)/(p2r.z - p1r.z);
				if (tmin >= 0.0 && tmin < t) {t = tmin; hit = 1;}
			}
			else {
				float const radius(0.5*(occlusion_only ? min(csz.x, csz.y) : max(csz.x, csz.y))); // use conservative radius unless this is an occlusion query
				point const cp1(cc.x, cc.y, i->z1()), cp2(cc.x, cc.y, i->z2());
				if (!line_int_cylinder(p1r, p2r, cp1, cp2, radius, radius, 1, tmin) || tmin > t) continue; // conservative for non-occlusion rays

				if (!occlusion_only && csz.x != csz.y) { // ellipse
					vector3d const delta(p2r - p1r);
					float const rx_inv_sq(1.0/(0.25*csz.x*csz.x)), ry_inv_sq(1.0/(0.25*csz.y*csz.y));
					float t_step(0.1*max(csz.x, csz.y)/delta.mag());

					for (unsigned n = 0; n < 10; ++n) { // use an interative approach
						if (point_in_ellipse_risq((p1r + tmin*delta), cc, rx_inv_sq, ry_inv_sq)) {hit = 1; tmin -= t_step;} else {tmin += t_step;}
						if (hit) {t_step *= 0.5;} // converge on hit point
					}
					if (!hit) continue; // not actually a hit
				} // end ellipse case
				t = tmin; hit = 1;
			}
		}
		else if (num_sides != 4) {
			building_draw_utils::calc_poly_pts(*this, *i, points);
			float const tz((i->z2() - p1r.z)/(p2r.z - p1r.z)); // t value at zval = top of cube

			if (tz >= 0.0 && tz < t) {
				float const xval(p1r.x + tz*(p2r.x - p1r.x)), yval(p1r.y + tz*(p2r.y - p1r.y));
				if (point_in_polygon_2d(xval, yval, &points.front(), points.size(), 0, 1)) {t = tz; hit = 1;} // XY plane test for vertical lines and top surface
			}
			if (!vert) { // test building sides
				point quad_pts[4]; // quads

				for (unsigned S = 0; S < num_sides; ++S) { // generate vertex data quads
					for (unsigned d = 0, ix = 0; d < 2; ++d) {
						point const &p(points[(S+d)%num_sides]);
						for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, i->d[2][d^e]);}
					}
					if (line_poly_intersect(p1r, p2r, quad_pts, 4, get_poly_norm(quad_pts), tmin) && tmin < t) {t = tmin; hit = 1;} // Note: untested
				} // for S
			}
		}
		else if (get_line_clip(p1r, p2r, i->d, tmin, tmax) && tmin < t) {t = tmin; hit = 1;} // cube

		if (hit) {
			if (occlusion_only) return 1; // early exit
			float const zval(p1.z + t*(p2.z - p1.z));
			coll = ((fabs(zval - i->d[2][1]) < 0.0001*i->get_dz()) ? 2 : 1); // test if clipped zval is close to the roof zval
			if (ret_any_pt) return coll;
		}
	} // for i
	if (occlusion_only) return 0;

	for (auto i = details.begin(); i != details.end(); ++i) {
		if (get_line_clip(p1r, p2r, i->d, tmin, tmax) && tmin < t) {t = tmin; coll = 3;} // details cube
	}
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
		if (line_poly_intersect(p1r, p2r, i->pts, i->npts, i->get_norm(), tmin) && tmin < t) {t = tmin; coll = 2;} // roof quad
	}
	return coll; // Note: no collisions with windows or doors, since they're colinear with walls
}

// Note: if xy_radius == 0.0, this is a point test; otherwise, it's an approximate vertical cylinder test
bool building_t::check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const {

	if (xy_radius == 0.0 && !bcube.contains_pt(pos)) return 0; // no intersection
	point pr(pos);
	if (is_rotated()) {do_xy_rotate(-rot_sin, rot_cos, bcube.get_cube_center(), pr);} // inverse rotate - negate the sine term
	
	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (pr.z > i->z2() || pr.z < i->z1()) continue; // no overlap in z

		if (use_cylinder_coll()) { // vertical cylinder
			point const cc(i->get_cube_center());
			vector3d const csz(i->get_size());
			float const dx(cc.x - pr.x), dy(cc.y - pr.y), rx(0.5*csz.x + xy_radius), ry(0.5*csz.y + xy_radius);
			if (dx*dx/(rx*rx) + dy*dy/(ry*ry) > 1.0) continue; // no intersection (below test should return true as well)
			return 1;
		}
		else if (num_sides != 4) {
			building_draw_utils::calc_poly_pts(*this, *i, points);

			if (xy_radius > 0.0) { // cylinder case: expand polygon by xy_radius; assumes a convex polygon
				point const center(i->get_cube_center());

				for (auto i = points.begin(); i != points.end(); ++i) {
					vector3d dir(*i - center);
					dir.z = 0.0; // only want XY component
					*i += dir*(xy_radius/dir.mag());
				}
			}
			if (point_in_polygon_2d(pr.x, pr.y, &points.front(), points.size(), 0, 1)) return 1; // XY plane test for top surface
		}
		else { // cube
			if (xy_radius > 0.0) {
				cube_t cube(*i);
				cube.expand_by(xy_radius);
				if (cube.contains_pt(pr)) return 1;
			}
			else if (i->contains_pt(pr)) return 1;
		}
	} // for i
	return 0;
}

void building_t::calc_bcube_from_parts() {
	assert(!parts.empty());
	bcube = parts[0];
	for (auto i = parts.begin()+1; i != parts.end(); ++i) {bcube.union_with_cube(*i);} // update bcube
}

void building_t::gen_geometry(int rseed1, int rseed2) {

	if (!is_valid()) return; // invalid building
	cube_t const base(parts.empty() ? bcube : parts.back());
	assert(base.is_strictly_normalized());
	parts.clear();
	details.clear();
	roof_tquads.clear();
	doors.clear();
	building_mat_t const &mat(get_material());
	rand_gen_t rgen;
	rgen.set_state(123+rseed1, 345*rseed2);
	ao_bcz2 = bcube.z2(); // capture z2 before union with roof and detail geometry (which increases building height)
	if (is_house) {gen_house(base, rgen); return;}

	// determine building shape (cube, cylinder, other)
	if (rgen.rand_probability(mat.round_prob)) {num_sides = MAX_CYLIN_SIDES;} // max number of sides for drawing rounded (cylinder) buildings
	else if (rgen.rand_probability(mat.cube_prob)) {num_sides = 4;} // cube
	else { // N-gon
		num_sides = mat.min_sides;
		if (mat.min_sides != mat.max_sides) {num_sides += (rgen.rand() % (1 + abs((int)mat.max_sides - (int)mat.min_sides)));}
	}
	bool const was_cube(is_cube()); // before num_sides increase due to ASF

	if (num_sides >= 6 && mat.max_fsa > 0.0) { // at least 6 sides
		flat_side_amt = max(0.0f, min(0.45f, rgen.rand_uniform(mat.min_fsa, mat.max_fsa)));
		if (flat_side_amt > 0.0 && rot_sin == 0.0) {start_angle = rgen.rand_uniform(0.0, TWO_PI);} // flat side, not rotated: add random start angle to break up uniformity
	}
	if ((num_sides == 3 || num_sides == 4 || num_sides == 6) && mat.max_asf > 0.0 && rgen.rand_probability(mat.asf_prob)) { // triangles/cubes/hexagons
		alt_step_factor = max(0.0f, min(0.99f, rgen.rand_uniform(mat.min_asf, mat.max_asf)));
		if (alt_step_factor > 0.0 && !(num_sides&1)) {half_offset = 1;} // chamfered cube/hexagon
		if (alt_step_factor > 0.0) {num_sides *= 2;}
	}

	// determine the number of levels and splits
	unsigned num_levels(mat.min_levels);
	
	if (mat.min_levels < mat.max_levels) { // have a range of levels
		if (was_cube || rgen.rand_bool()) {num_levels += rgen.rand() % (mat.max_levels - mat.min_levels + 1);} // only half of non-cubes are multilevel (unless min_level > 1)
	}
	if (mat.min_level_height > 0.0) {num_levels = max(mat.min_levels, min(num_levels, unsigned(bcube.get_size().z/mat.min_level_height)));}
	num_levels = max(num_levels, 1U); // min_levels can be zero to apply more weight to 1 level buildings
	bool const do_split(num_levels < 4 && is_cube() && rgen.rand_probability(mat.split_prob)); // don't split buildings with 4 or more levels, or non-cubes

	if (num_levels == 1) { // single level
		if (do_split) {split_in_xy(base, rgen);} // generate L, T, or U shape
		else { // single part, entire cube/cylinder
			parts.push_back(base);
			if ((rgen.rand()&3) != 0) {gen_sloped_roof(rgen);} // 75% chance
			gen_details(rgen);
		}
		return; // for now the bounding cube
	}
	// generate building levels and splits
	parts.resize(num_levels);
	float const height(base.get_dz()), dz(height/num_levels);
	assert(height > 0.0);

	if (!do_split && (rgen.rand()&3) < (was_cube ? 2 : 3)) { // oddly shaped multi-sided overlapping sections (50% chance for cube buildings and 75% chance for others)
		point const llc(base.get_llc()), sz(base.get_size());

		for (unsigned i = 0; i < num_levels; ++i) { // generate overlapping cube levels
			cube_t &bc(parts[i]);
			bc.z1() = base.z1(); // z1
			bc.z2() = base.z1() + (i+1)*dz; // z2
			if (i > 0) {bc.z2() += dz*rgen.rand_uniform(-0.5, 0.5); bc.z2() = min(bc.z2(), base.z2());}
			float const min_edge_mode(mat.no_city ? 0.04*i : 0.0); // prevent z-fighting on non-city building windows (stretched texture)

			for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to generate a cube that doesn't contain any existing cubes (can occasionally still fail)
				for (unsigned d = 0; d < 2; ++d) { // x,y
					bc.d[d][0] = base.d[d][0] + max(rgen.rand_uniform(-0.2, 0.45), min_edge_mode)*sz[d];
					bc.d[d][1] = base.d[d][1] - max(rgen.rand_uniform(-0.2, 0.45), min_edge_mode)*sz[d];
				}
				assert(bc.is_strictly_normalized());
				bool contains(0);
				for (unsigned j = 0; j < i; ++j) {contains |= bc.contains_cube(parts[j]);}
				if (!contains) break; // success
			} // for n
		} // for i
		calc_bcube_from_parts(); // update bcube
		gen_details(rgen);
		return;
	}
	for (unsigned i = 0; i < num_levels; ++i) {
		cube_t &bc(parts[i]);
		if (i == 0) {bc = base;} // use full building footprint
		else {
			cube_t const &prev(parts[i-1]);
			float const shift_mult(was_cube ? 1.0 : 0.5); // half the shift for non-cube buildings

			for (unsigned d = 0; d < 2; ++d) {
				float const len(prev.d[d][1] - prev.d[d][0]), bc_len(bcube.d[d][1] - bcube.d[d][0]);
				bool const inv(rgen.rand_bool());

				for (unsigned E = 0; E < 2; ++E) {
					bool const e((E != 0) ^ inv); // no dir favoritism for 20% check
					float delta(0.0);
					if (rgen.rand()&3) {delta = shift_mult*rgen.rand_uniform(0.1, 0.4);} // 25% chance of no shift, 75% chance of 20-40% shift
					bc.d[d][e] = prev.d[d][e] + (e ? -delta : delta)*len;
					if (bc.d[d][1] - bc.d[d][0] < (0.2/shift_mult)*bc_len) {bc.d[d][e] = prev.d[d][e];} // if smaller than 20% base width, revert the change
				}
			}
			bc.z1() = prev.z2(); // z1
		}
		bc.z2() = bc.z1() + dz; // z2
		bc.normalize(); // handle XY inversion due to shift
	} // for i
	for (unsigned i = 1; i < num_levels; ++i) {
		float const ddz(rgen.rand_uniform(-0.35*dz, 0.35*dz)); // random shift in z height
		parts[i  ].z1() += ddz;
		parts[i-1].z2() += ddz;
	}
	if (do_split) { // generate L, T, or U shape
		cube_t const split_cube(parts.back());
		parts.pop_back();
		split_in_xy(split_cube, rgen);
	}
	else {
		if ((rgen.rand()&3) != 0) {gen_sloped_roof(rgen);} // 67% chance
		if (num_levels <= 3) {gen_details(rgen);}
	}
}

bool get_largest_xy_dim(cube_t const &c) {return (c.dy() > c.dx());}

void building_t::gen_house(cube_t const &base, rand_gen_t &rgen) {

	assert(parts.empty());
	int const type(rgen.rand()%3); // 0=single cube, 1=L-shape, 2=two-part
	unsigned force_dim[2] = {2}; // force roof dim to this value, per part; 2 = unforced/auto
	bool skip_last_roof(0);
	num_sides = 4;
	if (type != 0) {parts.reserve(4);} // two house sections + porch roof + porch support (upper bound)
	parts.push_back(base);
	// add a door
	bool const gen_door(global_building_params.windows_enabled());
	float door_height(0.45/(get_material().wind_yscale*global_building_params.get_window_ty())); // set height based on window spacing
	float door_center(0.0), door_pos(0.0);
	door_dim = rgen.rand_bool();

	if (type != 0) { // multi-part house
		parts.push_back(base); // add second part
		bool const dir(rgen.rand_bool()); // in dim
		float const split(rgen.rand_uniform(0.4, 0.6)*(dir  ? -1.0 : 1.0));
		float delta_height(0.0), shrink[2] = {0.0};
		bool dim(0), dir2(0);

		if (type == 1) { // L-shape
			dir2         = rgen.rand_bool(); // in !dim
			dim          = rgen.rand_bool();
			shrink[dir2] = rgen.rand_uniform(0.4, 0.6)*(dir2 ? -1.0 : 1.0);
			delta_height = max(0.0f, rand_uniform(-0.1, 0.5));
		}
		else if (type == 2) { // two-part
			dim          = get_largest_xy_dim(base); // choose longest dim
			delta_height = rand_uniform(0.1, 0.5);
			
			for (unsigned d = 0; d < 2; ++d) {
				if (rgen.rand_bool()) {shrink[d] = rgen.rand_uniform(0.2, 0.35)*(d ? -1.0 : 1.0);}
			}
		}
		vector3d const sz(base.get_size());
		parts[0].d[ dim][ dir] += split*sz[dim]; // split in dim
		parts[1].d[ dim][!dir]  = parts[0].d[dim][dir];
		cube_t const pre_shrunk_p1(parts[1]); // save for use in details below
		for (unsigned d = 0; d < 2; ++d) {parts[1].d[!dim][d] += shrink[d]*sz[!dim];} // shrink this part in the other dim
		parts[1].z2() -= delta_height*parts[1].dz(); // lower height
		if (type == 1 && rgen.rand_bool()) {force_dim[0] = !dim; force_dim[1] = dim;} // L-shape, half the time
		else if (type == 2) {force_dim[0] = force_dim[1] = dim;} // two-part - force both parts to have roof along split dim
		int const detail_type((type == 1) ? (rgen.rand()%3) : 0); // 0=none, 1=porch, 2=detatched garage/shed
		door_dir = ((door_dim == dim) ? dir : dir2); // if we have a porch/shed/garage, put the door on that side
		if (door_dim == dim && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way

		if (detail_type != 0) { // add details to L-shaped house
			cube_t c(pre_shrunk_p1);
			c.d[!dim][!dir2] = parts[1].d[!dim][dir2]; // other half of the shrunk part1
			float const dist1((c.d[!dim][!dir2] - base.d[!dim][dir2])*rgen.rand_uniform(0.4, 0.6));
			float const dist2((c.d[ dim][!dir ] - base.d[ dim][dir ])*rgen.rand_uniform(0.4, 0.6));
			float const height(rgen.rand_uniform(0.55, 0.7)*parts[1].dz());
			
			if (gen_door) { // add door in interior of L, centered under porch roof (if it exists, otherwise where it would be)
				door_center = 0.5*(c.d[!door_dim][0] + c.d[!door_dim][1] + ((door_dim == dim) ? dist1 : dist2));
				door_pos    = c.d[door_dim][!door_dir];
				door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
				min_eq(door_height, 0.95f*height);
			}
			if (detail_type == 1) { // porch
				float const width(0.05*(fabs(dist1) + fabs(dist2))); // width of support pillar
				c.d[!dim][dir2] += dist1; // move away from bcube edge
				c.d[ dim][ dir] += dist2; // move away from bcube edge
				c.z1() += height; // move up
				c.z2()  = c.z1() + 0.05*parts[1].dz();
				parts.push_back(c); // porch roof
				c.z2() = c.z1();
				c.z1() = pre_shrunk_p1.z1(); // support pillar
				c.d[!dim][!dir2] = c.d[!dim][dir2] + (dir2 ? -1.0 : 1.0)*width;
				c.d[ dim][!dir ] = c.d[ dim][ dir] + (dir  ? -1.0 : 1.0)*width;
				skip_last_roof = 1;
			}
			else if (detail_type == 2) { // detatched garage/shed
				c.d[!dim][dir2 ]  = base.d[!dim][dir2]; // shove it into the opposite corner of the bcube
				c.d[ dim][dir  ]  = base.d[ dim][dir ]; // shove it into the opposite corner of the bcube
				c.d[!dim][!dir2] -= dist1; // move away from bcube edge
				c.d[ dim][!dir ] -= dist2; // move away from bcube edge
				c.z2() = c.z1() + min(min(c.dx(), c.dy()), height); // no taller than x or y size; Note: z1 same as part1
			}
			parts.push_back(c);
		} // end house details
		calc_bcube_from_parts(); // maybe calculate a tighter bounding cube
	} // end type != 0
	else if (gen_door) { // single cube house
		door_dir  = rgen.rand_bool(); // select a random dir
		door_part = 0; // only one part
	}
	if (gen_door) { // add door
		if (door_center == 0.0) { // door not yet calculated
			cube_t const &c(parts[0]); // add door to first part of house
			float const offset(rgen.rand_uniform(0.25, 0.75));
			door_center = offset*c.d[!door_dim][0] + (1.0 - offset)*c.d[!door_dim][1];
			door_pos    = c.d[door_dim][door_dir];
			door_part   = 0;
		}
		float const door_half_width(0.25*door_height);
		cube_t door;
		door.z1() = base.z1(); // same bottom as house
		door.z2() = door.z1() + door_height;
		door.d[ door_dim][!door_dir] = door_pos + 0.005*base.dz()*(door_dir ? 1.0 : -1.0); // move slightly away from the house to prevent z-fighting
		door.d[ door_dim][ door_dir] = door.d[door_dim][!door_dir]; // make zero size in this dim
		door.d[!door_dim][0] = door_center - door_half_width; // left
		door.d[!door_dim][1] = door_center + door_half_width; // right
		add_door(door, door_dim, door_dir);
	}
	float const peak_height(rgen.rand_uniform(0.15, 0.5)); // same for all parts

	for (auto i = parts.begin(); (i + skip_last_roof) != parts.end(); ++i) {
		unsigned const fdim(force_dim[i - parts.begin()]);
		bool const dim((fdim < 2) ? fdim : get_largest_xy_dim(*i)); // use longest side if not forced
		gen_peaked_roof(*i, peak_height, dim);
	}
	add_roof_to_bcube();
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for roof
}

void building_t::add_door(cube_t const &c, bool dim, bool dir) {

	vector3d const sz(c.get_size());
	assert(sz[dim] == 0.0 && sz[!dim] > 0.0 && sz.z > 0.0);
	tquad_with_ix_t door(4, tquad_with_ix_t::TYPE_DOOR); // quad
	door.pts[0].z = door.pts[1].z = c.z1(); // bottom
	door.pts[2].z = door.pts[3].z = c.z2(); // top
	door.pts[0][!dim] = door.pts[3][!dim] = c.d[!dim][ dir]; //  dir side
	door.pts[1][!dim] = door.pts[2][!dim] = c.d[!dim][!dir]; // !dir side
	door.pts[0][ dim] = door.pts[1][ dim] = door.pts[2][ dim] = door.pts[3][ dim] = c.d[dim][0] + 0.01*sz[!dim]*(dir ? 1.0 : -1.0); // move away from wall slightly
	doors.push_back(door);
}

void building_t::gen_peaked_roof(cube_t const &top, float peak_height, bool dim) { // roof made from two sloped quads

	float const width(dim ? top.get_dx() : top.get_dy()), roof_dz(min(peak_height*width, top.get_dz()));
	float const z1(top.z2()), z2(z1 + roof_dz), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	point pts[6] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(x1, y1, z2), point(x2, y2, z2)};
	if (dim == 0) {pts[4].y = pts[5].y = 0.5*(y1 + y2);} // yc
	else          {pts[4].x = pts[5].x = 0.5*(x1 + x2);} // xc
	unsigned const qixs[2][2][4] = {{{0,3,5,4}, {4,5,2,1}}, {{0,4,5,1}, {4,3,2,5}}}; // 2 quads
	roof_tquads.reserve(4); // 2 roof quads + 2 side triangles

	// TODO: extend outside the wall a small amount? may require updating bcube for drawing
	for (unsigned n = 0; n < 2; ++n) { // roof
		tquad_t tquad(4); // quad
		UNROLL_4X(tquad.pts[i_] = pts[qixs[dim][n][i_]];)
		roof_tquads.emplace_back(tquad, tquad_with_ix_t::TYPE_ROOF); // tag as roof
	}
	unsigned const tixs[2][2][3] = {{{1,0,4}, {3,2,5}}, {{0,3,4}, {2,1,5}}}; // 2 triangles

	for (unsigned n = 0; n < 2; ++n) { // triangle section/wall from z1 up to roof
		tquad_t tquad(3); // triangle
		UNROLL_3X(tquad.pts[i_] = pts[tixs[dim][n][i_]];)
		roof_tquads.emplace_back(tquad, tquad_with_ix_t::TYPE_WALL); // tag as wall
	}
}

void building_t::gen_details(rand_gen_t &rgen) { // for the roof

	unsigned const num_blocks(roof_tquads.empty() ? (rgen.rand() % 9) : 0); // 0-8; 0 if there are roof quads (houses, etc.)
	has_antenna = (rgen.rand() & 1);
	details.resize(num_blocks + has_antenna);
	assert(!parts.empty());
	if (details.empty()) return; // nothing to do
	cube_t const &top(parts.back()); // top/last part

	if (num_blocks > 0) {
		float const xy_sz(top.get_size().xy_mag());
		float const height_scale(0.0035*(top.dz() + bcube.dz())); // based on avg height of current section and entire building
		cube_t rbc(top);
		vector<point> points; // reused across calls
		
		for (unsigned i = 0; i < num_blocks; ++i) {
			cube_t &c(details[i]);
			float const height(height_scale*rgen.rand_uniform(1.0, 4.0));

			while (1) {
				c.set_from_point(point(rgen.rand_uniform(rbc.x1(), rbc.x2()), rgen.rand_uniform(rbc.y1(), rbc.y2()), 0.0));
				c.expand_by(vector3d(xy_sz*rgen.rand_uniform(0.01, 0.08), xy_sz*rgen.rand_uniform(0.01, 0.06), 0.0));
				if (!rbc.contains_cube_xy(c)) continue; // not contained
				if (is_simple_cube()) break; // success/done
				bool contained(1);

				for (unsigned i = 0; i < 4; ++i) { // check cylinder/ellipse
					point const pt(c.d[0][i&1], c.d[1][i>>1], 0.0); // XY only
					if (!check_part_contains_pt_xy(rbc, pt, points)) {contained = 0; break;}
				}
				if (contained) break; // success/done
			} // end while
			c.z1() = top.z2(); // z1
			c.z2() = top.z2() + height; // z2
		} // for i
	}
	if (has_antenna) { // add antenna
		float const radius(0.003*rgen.rand_uniform(1.0, 2.0)*(top.get_dx() + top.get_dy()));
		float const height(rgen.rand_uniform(0.25, 0.5)*top.get_dz());
		cube_t &antenna(details.back());
		antenna.set_from_point(top.get_cube_center());
		antenna.expand_by(vector3d(radius, radius, 0.0));
		antenna.z1() = top.z2(); // z1
		antenna.z2() = bcube.z2() + height; // z2 (use bcube to include sloped roof)
	}
	for (auto i = details.begin(); i != details.end(); ++i) {max_eq(bcube.z2(), i->z2());} // extend bcube z2 to contain details
	if (roof_tquads.empty()) {gen_grayscale_detail_color(rgen, 0.2, 0.6);} // for antenna and roof
}

void building_t::gen_sloped_roof(rand_gen_t &rgen) { // Note: currently not supported for rotated buildings

	assert(!parts.empty());
	if (!is_simple_cube()) return; // only simple cubes are handled
	cube_t const &top(parts.back()); // top/last part
	float const peak_height(rgen.rand_uniform(0.2, 0.5));
	float const wmin(min(top.get_dx(), top.get_dy())), z1(top.z2()), z2(z1 + peak_height*wmin), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	point const pts[5] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(0.5*(x1 + x2), 0.5*(y1 + y2), z2)};
	float const d1(rgen.rand_uniform(0.0, 0.8));

	if (d1 < 0.2) { // pointed roof with 4 sloped triangles
		unsigned const ixs[4][3] = {{1,0,4}, {3,2,4}, {0,3,4}, {2,1,4}};
		roof_tquads.resize(4);
	
		for (unsigned n = 0; n < 4; ++n) {
			roof_tquads[n].npts = 3; // triangles
			UNROLL_3X(roof_tquads[n].pts[i_] = pts[ixs[n][i_]];)
		}
	}
	else { // flat roof with center quad and 4 surrounding sloped quads
		point const center((1.0 - d1)*pts[4]);
		point pts2[8];
		for (unsigned n = 0; n < 4; ++n) {pts2[n] = pts[n]; pts2[n+4] = d1*pts[n] + center;}
		unsigned const ixs[5][4] = {{4,7,6,5}, {0,4,5,1}, {3,2,6,7}, {0,3,7,4}, {2,1,5,6}}; // add the flat quad first, which works better for sphere intersections
		roof_tquads.resize(5);

		for (unsigned n = 0; n < 5; ++n) {
			roof_tquads[n].npts = 4; // quads
			UNROLL_4X(roof_tquads[n].pts[i_] = pts2[ixs[n][i_]];)
		}
	}
	add_roof_to_bcube();
	//max_eq(bcube.z2(), z2);
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for antenna and roof
}

void building_t::add_roof_to_bcube() {
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {i->update_bcube(bcube);} // technically should only need to update z2
}
void building_t::gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax) {
	float const cscale(rgen.rand_uniform(imin, imax));
	detail_color = colorRGBA(cscale, cscale, cscale, 1.0);
}

bool check_tile_smap(bool shadow_only) {
	return (!shadow_only && world_mode == WMODE_INF_TERRAIN && shadow_map_enabled());
}

void building_t::get_all_drawn_verts(building_draw_t &bdraw) const {

	if (!is_valid()) return; // invalid building
	building_mat_t const &mat(get_material());

	for (auto i = parts.begin(); i != parts.end(); ++i) { // multiple cubes/parts/levels
		bdraw.add_section(*this, *i, bcube, (is_house ? i->z2() : ao_bcz2), mat.side_tex, side_color, 3, 0, 0, 0, 0); // XY
		bool const skip_top(!roof_tquads.empty() && (is_house || i+1 == parts.end())); // don't add the flat roof for the top part in this case
		bool const is_stacked(!is_house && num_sides == 4 && i->d[2][0] > bcube.d[2][0]); // skip the bottom of stacked cubes
		if (is_stacked && skip_top) continue; // no top/bottom to draw
		bdraw.add_section(*this, *i, bcube, (is_house ? i->z2() : ao_bcz2), mat.roof_tex, roof_color, 4, is_stacked, skip_top, 0, 0); // only Z dim
	}
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
		bool const is_wall(i->type == tquad_with_ix_t::TYPE_WALL);
		bdraw.add_tquad(*this, *i, bcube, (is_wall ? mat.side_tex : mat.roof_tex), (is_wall ? side_color : roof_color)); // use type to select roof vs. side texture
	}
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		bdraw.add_tquad(*this, *i, bcube, tid_nm_pair_t(building_window_gen.get_door_tid(), -1, 1.0, 1.0), WHITE);
	}
	for (auto i = details.begin(); i != details.end(); ++i) { // draw roof details
		building_geom_t bg(4, rot_sin, rot_cos); // cube
		bg.is_pointed = (has_antenna && i+1 == details.end()); // draw antenna as a point
		bdraw.add_section(bg, *i, bcube, ao_bcz2, mat.roof_tex.get_scaled_version(0.5), detail_color*(bg.is_pointed ? 0.5 : 1.0), 7, 1, 0, 1, 0); // all dims, skip_bottom, no AO
	}
}

void building_t::get_all_drawn_window_verts(building_draw_t &bdraw, bool lights_pass) const {

	if (!is_valid() || !global_building_params.windows_enabled()) return; // invalid building or no windows
	building_mat_t const &mat(get_material());
	if (lights_pass ? !mat.add_wind_lights : !mat.add_windows) return; // no windows for this material
	int const window_tid(building_window_gen.get_window_tid());
	if (window_tid < 0) return; // not allocated - error?
	if (mat.wind_xscale == 0.0 || mat.wind_yscale == 0.0) return; // no windows for this material?
	tid_nm_pair_t tex(window_tid, -1, mat.wind_xscale*global_building_params.get_window_tx(), mat.wind_yscale*global_building_params.get_window_ty(), mat.wind_xoff, mat.wind_yoff);
	colorRGBA color;

	if (lights_pass) { // slight yellow-blue tinting using bcube x1/y1 as a hash
		float const tint(0.2*fract(100.0*(bcube.x1() + bcube.y1())));
		color = colorRGBA((1.0 - tint), (1.0 - tint), (0.8 + tint), 1.0);
	}
	else {color = mat.window_color;}
	bool const clip_windows(mat.no_city); // only clip non-city windows; city building windows tend to be aligned with the building textures (maybe should be a material option?)

	for (auto i = parts.begin(); i != parts.end(); ++i) { // multiple cubes/parts/levels
		float const door_ztop((!doors.empty() && (i - parts.begin()) == unsigned(door_part)) ? doors.front().pts[2].z : 0.0);
		bdraw.add_section(*this, *i, bcube, ao_bcz2, tex, color, 3, 0, 0, 1, clip_windows, door_ztop); // XY, no_ao=1
	}
}


class building_creator_t {

	unsigned grid_sz;
	vector3d range_sz, range_sz_inv, max_extent;
	cube_t range, buildings_bcube;
	rand_gen_t rgen;
	vector<building_t> buildings;
	vector<vector<unsigned>> bix_by_plot; // cached for use with pedestrian collisions
	building_draw_t building_draw, building_draw_vbo, building_draw_windows, building_draw_wind_lights;
	vector<point> points; // reused temporary

	struct grid_elem_t {
		vector<unsigned> ixs;
		cube_t bcube;
		grid_elem_t() : bcube(all_zeros) {}
		void add(cube_t const &c, unsigned ix) {
			if (ixs.empty()) {bcube = c;} else {bcube.union_with_cube(c);}
			ixs.push_back(ix);
		}
	};
	vector<grid_elem_t> grid, grid_by_tile;

	grid_elem_t &get_grid_elem(unsigned gx, unsigned gy) {
		assert(gx < grid_sz && gy < grid_sz);
		return grid[gy*grid_sz + gx];
	}
	grid_elem_t const &get_grid_elem(unsigned gx, unsigned gy) const {
		assert(gx < grid_sz && gy < grid_sz);
		return grid[gy*grid_sz + gx];
	}
	struct bix_by_x1 {
		vector<building_t> const &buildings;
		bix_by_x1(vector<building_t> const &buildings_) : buildings(buildings_) {}
		bool operator()(unsigned const a, unsigned const b) const {return (buildings[a].bcube.x1() < buildings[b].bcube.x1());}
	};
	unsigned get_grid_ix(point pos) const {
		range.clamp_pt_xy(pos);
		unsigned gxy[2];
		for (unsigned d = 0; d < 2; ++d) {
			float const v((pos[d] - range.d[d][0])*range_sz_inv[d]);
			gxy[d] = unsigned(v*(grid_sz-1));
			assert(gxy[d] < grid_sz);
		}
		return (gxy[1]*grid_sz + gxy[0]);
	}
	void get_grid_range(cube_t const &bcube, unsigned ixr[2][2]) const { // {lo,hi}x{x,y}
		point llc(bcube.get_llc()), urc(bcube.get_urc());
		range.clamp_pt_xy(llc);
		range.clamp_pt_xy(urc);
		for (unsigned d = 0; d < 2; ++d) {
			float const v1((llc[d] - range.d[d][0])*range_sz_inv[d]), v2((urc[d] - range.d[d][0])*range_sz_inv[d]);
			ixr[0][d] = unsigned(v1*(grid_sz-1));
			ixr[1][d] = unsigned(v2*(grid_sz-1));
			assert(ixr[0][d] < grid_sz && ixr[1][d] < grid_sz);
		}
	}
	void add_to_grid(cube_t const &bcube, unsigned bix) {
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				get_grid_elem(x, y).add(bcube, bix);
			}
		}
	}
	bool check_for_overlaps(vector<unsigned> const &ixs, cube_t const &test_bc, building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const {
		for (auto i = ixs.begin(); i != ixs.end(); ++i) {
			building_t const &ob(get_building(*i));
			if (test_bc.intersects_xy(ob.bcube) && ob.check_bcube_overlap_xy(b, expand_rel, expand_abs, points)) return 1;
		}
		return 0;
	}

	void build_grid_by_tile(bool single_tile) {
		grid_by_tile.clear();

		if (single_tile || world_mode != WMODE_INF_TERRAIN) { // not used in this mode - add all buildings to the first tile
			grid_by_tile.resize(1);
			grid_by_tile.front().ixs.reserve(buildings.size());
			for(unsigned bix = 0; bix < buildings.size(); ++bix) {grid_by_tile.front().add(buildings[bix].bcube, bix);}
			return;
		}
		//timer_t timer("build_grid_by_tile");
		map<uint64_t, unsigned> tile_to_gbt;

		for(unsigned bix = 0; bix < buildings.size(); ++bix) {
			unsigned gix;
			cube_t const &bcube(buildings[bix].bcube);
			uint64_t const tile_id(get_tile_id_containing_point_no_xyoff(bcube.get_cube_center()));
			auto it(tile_to_gbt.find(tile_id));

			if (it == tile_to_gbt.end()) { // new element
				gix = grid_by_tile.size();
				grid_by_tile.push_back(grid_elem_t());
				tile_to_gbt[tile_id] = gix;
			}
			else { // existing element
				gix = it->second;
				assert(gix < grid_by_tile.size());
			}
			grid_by_tile[gix].add(bcube, bix);
		} // for bix
	}

	bool check_valid_building_placement(building_params_t const &params, building_t const &b, vect_cube_t const &avoid_bcubes,
		float min_building_spacing, unsigned plot_ix, bool non_city_only, bool use_city_plots, bool check_plot_coll)
	{
		float const expand_val(b.is_rotated() ? 0.05 : 0.1); // expand by 5-10% (relative - multiplied by building size)
		vector3d const b_sz(b.bcube.get_size());
		vector3d expand(expand_val*b_sz);
		for (unsigned d = 0; d < 2; ++d) {max_eq(expand[d], min_building_spacing);} // ensure the min building spacing (only applies to the current building)
		cube_t test_bc(b.bcube);
		test_bc.expand_by_xy(expand);

		if (use_city_plots) {
			assert(plot_ix < bix_by_plot.size());
			if (check_for_overlaps(bix_by_plot[plot_ix], test_bc, b, expand_val, min_building_spacing, points)) return 0;
			bix_by_plot[plot_ix].push_back(buildings.size());
		}
		else if (check_plot_coll && has_bcube_int_xy(test_bc, avoid_bcubes, params.sec_extra_spacing)) { // extra expand val
			return 0;
		}
		else {
			float const extra_spacing(non_city_only ? params.sec_extra_spacing : 0.0); // absolute value of expand
			test_bc.expand_by_xy(extra_spacing);
			unsigned ixr[2][2];
			get_grid_range(test_bc, ixr);

			for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
				for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
					grid_elem_t const &ge(get_grid_elem(x, y));
					if (!test_bc.intersects_xy(ge.bcube)) continue;
					if (check_for_overlaps(ge.ixs, test_bc, b, expand_val, max(min_building_spacing, extra_spacing), points)) {return 0;}
				} // for x
			} // for y
		}
		return 1;
	}

	struct building_cand_t : public building_t {
		vect_cube_t &temp_parts;
		building_cand_t(vect_cube_t &temp_parts_) : temp_parts(temp_parts_) {temp_parts.clear(); parts.swap(temp_parts);} // parts takes temp_parts memory
		~building_cand_t() {parts.swap(temp_parts);} // memory returned from parts to temp_parts
	};

public:
	building_creator_t() : grid_sz(1), max_extent(zero_vector) {}
	bool empty() const {return buildings.empty();}
	void clear() {buildings.clear(); grid.clear();}
	vector3d const &get_max_extent() const {return max_extent;}
	building_t const &get_building(unsigned ix) const {assert(ix < buildings.size()); return buildings[ix];}
	cube_t const &get_building_bcube(unsigned ix) const {return get_building(ix).bcube;}
	cube_t const &get_bcube() const {return buildings_bcube;}
	bool is_visible(vector3d const &xlate) const {return (!empty() && camera_pdu.cube_visible(buildings_bcube + xlate));}
	
	bool get_building_hit_color(point const &p1, point const &p2, colorRGBA &color) const {
		float t(0.0); // unused
		unsigned hit_bix(0);
		unsigned const ret(check_line_coll(p1, p2, t, hit_bix, 0)); // 0=no hit, 1=hit side, 2=hit roof
		if (ret == 0) return 0;
		building_t const &b(get_building(hit_bix));
		switch (ret) {
		case 1: color = b.get_avg_side_color  (); break;
		case 2: color = b.get_avg_roof_color  (); break;
		case 3: color = b.get_avg_detail_color(); break;
		default: assert(0);
		}
		return 1;
	}
	void gen(building_params_t const &params, bool city_only, bool non_city_only, bool is_tile, int rseed=123) {
		assert(!(city_only && non_city_only));
		clear();
		if (params.tt_only && world_mode != WMODE_INF_TERRAIN) return;
		if (params.gen_inf_buildings() && !is_tile) return; // not added here
		vector<unsigned> const &mat_ix_list(params.get_mat_list(city_only, non_city_only));
		if (params.materials.empty() || mat_ix_list.empty()) return; // no materials
		timer_t timer("Gen Buildings");
		float const def_water_level(get_water_z_height()), min_building_spacing(get_min_obj_spacing());
		vector3d const offset(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0);
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? offset : zero_vector); // cancel out xoff2/yoff2 translate
		vector3d const delta_range((world_mode == WMODE_INF_TERRAIN) ? zero_vector : offset);
		range = params.materials[mat_ix_list.front()].pos_range; // range is union over all material ranges
		for (auto i = mat_ix_list.begin()+1; i != mat_ix_list.end(); ++i) {range.union_with_cube(params.materials[*i].pos_range);}
		range     += delta_range;
		range_sz   = range.get_size(); // Note: place_radius is relative to range cube center
		max_extent = zero_vector;
		assert(range_sz.x > 0.0 && range_sz.y > 0.0);
		UNROLL_2X(range_sz_inv[i_] = 1.0/range_sz[i_];) // xy only
		if (!is_tile) {buildings.reserve(params.num_place);}
		grid_sz = (is_tile ? 4 : 32); // tiles are small enough that they don't need grids
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_tries(0), num_gen(0), num_skip(0);
		if (rseed == 0) {rseed = 123;} // 0 is a bad value
		rgen.set_state(rand_gen_index, rseed); // update when mesh changes, otherwise determinstic
		vect_cube_with_zval_t city_plot_bcubes;
		vect_cube_t avoid_bcubes;
		if (city_only) {get_city_plot_bcubes(city_plot_bcubes);} // Note: assumes approx equal area for placement distribution
		
		if (non_city_only) {
			get_city_bcubes(avoid_bcubes);
			get_city_road_bcubes(avoid_bcubes, 1); // connector roads only
			get_all_model_bcubes(avoid_bcubes);
			expand_cubes_by_xy(avoid_bcubes, get_road_max_width());
		}
		bool const use_city_plots(!city_plot_bcubes.empty()), check_plot_coll(!avoid_bcubes.empty());
		bix_by_plot.resize(city_plot_bcubes.size());
		point center(all_zeros);
		unsigned num_consec_fail(0), max_consec_fail(0);
		vect_cube_t temp_parts;

		for (unsigned i = 0; i < params.num_place; ++i) {
			bool success(0);

			for (unsigned n = 0; n < params.num_tries; ++n) { // 10 tries to find a non-overlapping building placement
				building_cand_t b(temp_parts);
				b.mat_ix = params.choose_rand_mat(rgen, city_only, non_city_only); // set material
				building_mat_t const &mat(b.get_material());
				cube_t pos_range;
				unsigned plot_ix(0);
				
				if (use_city_plots) { // select a random plot, if available
					plot_ix   = rgen.rand()%city_plot_bcubes.size();
					pos_range = city_plot_bcubes[plot_ix];
					center.z  = city_plot_bcubes[plot_ix].zval; // optimization: take zval from plot rather than calling get_exact_zval()
					pos_range.expand_by_xy(-min_building_spacing); // force min spacing between building and edge of plot
				}
				else {
					pos_range = mat.pos_range + delta_range;
				}
				vector3d const pos_range_sz(pos_range.get_size());
				assert(pos_range_sz.x > 0.0 && pos_range_sz.y > 0.0);
				point const place_center(pos_range.get_cube_center());
				bool keep(0);
				++num_tries;

				for (unsigned m = 0; m < params.num_tries; ++m) {
					for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(pos_range.d[d][0], pos_range.d[d][1]);} // x,y
					if (is_tile || mat.place_radius == 0.0 || dist_xy_less_than(center, place_center, mat.place_radius)) {keep = 1; break;} // place_radius ignored for tiles
				}
				if (!keep) continue; // placement failed, skip
				b.is_house = (mat.house_prob > 0.0 && rgen.rand_float() < mat.house_prob);
				float const size_scale(b.is_house ? mat.gen_size_scale(rgen) : 1.0);
				
				for (unsigned d = 0; d < 2; ++d) { // x,y
					float const sz(0.5*size_scale*rgen.rand_uniform(min(mat.sz_range.d[d][0], 0.3f*pos_range_sz[d]),
						                                            min(mat.sz_range.d[d][1], 0.5f*pos_range_sz[d]))); // use pos range size for max
					b.bcube.d[d][0] = center[d] - sz;
					b.bcube.d[d][1] = center[d] + sz;
				}
				if ((use_city_plots || is_tile) && !pos_range.contains_cube_xy(b.bcube)) continue; // not completely contained in plot/tile (pre-rot)
				if (!use_city_plots) {b.gen_rotation(rgen);} // city plots are Manhattan (non-rotated) - must rotate before bcube checks below
				if (is_tile && !pos_range.contains_cube_xy(b.bcube)) continue; // not completely contained in tile
				if (start_in_inf_terrain && b.bcube.contains_pt_xy(get_camera_pos())) continue; // don't place a building over the player appearance spot
				if (!check_valid_building_placement(params, b, avoid_bcubes, min_building_spacing, plot_ix, non_city_only, use_city_plots, check_plot_coll)) continue; // check overlap
				++num_gen;
				if (!use_city_plots) {center.z = get_exact_zval(center.x+xlate.x, center.y+xlate.y);} // only calculate when needed
				float const z_sea_level(center.z - def_water_level);
				if (z_sea_level < 0.0) break; // skip underwater buildings, failed placement
				if (z_sea_level < mat.min_alt || z_sea_level > mat.max_alt) break; // skip bad altitude buildings, failed placement
				float const hmin(use_city_plots ? pos_range.z1() : 0.0), hmax(use_city_plots ? pos_range.z2() : 1.0);
				assert(hmin <= hmax);
				float const height_range(mat.sz_range.dz());
				assert(height_range >= 0.0);
				float const z_size_scale(size_scale*(b.is_house ? rgen.rand_uniform(0.6, 0.8) : 1.0)); // make houses slightly shorter on average to offset extra height added by roof
				float const height_val(z_size_scale*(mat.sz_range.z1() + height_range*rgen.rand_uniform(hmin, hmax)));
				assert(height_val > 0.0);
				b.set_z_range(center.z, (center.z + 0.5*height_val));
				assert(b.bcube.is_strictly_normalized());
				mat.side_color.gen_color(b.side_color, rgen);
				mat.roof_color.gen_color(b.roof_color, rgen);
				add_to_grid(b.bcube, buildings.size());
				vector3d const sz(b.bcube.get_size());
				float const mult[3] = {0.5, 0.5, 1.0}; // half in X,Y and full in Z
				UNROLL_3X(max_extent[i_] = max(max_extent[i_], mult[i_]*sz[i_]);)
				if (buildings.empty()) {buildings_bcube = b.bcube;} else {buildings_bcube.union_with_cube(b.bcube);}
				buildings.push_back(b);
				success = 1;
				break; // done
			} // for n
			if (success) {num_consec_fail = 0;}
			else {
				++num_consec_fail;
				max_eq(max_consec_fail, num_consec_fail);

				if (num_consec_fail >= (is_tile ? 100U : 5000U)) { // too many failures - give up
					if (!is_tile) {cout << "Failed to place a building after " << num_consec_fail << " tries, giving up after " << i << " iterations" << endl;}
					break;
				}
			}
		} // for i
		bix_by_x1 cmp_x1(buildings);
		for (auto i = bix_by_plot.begin(); i != bix_by_plot.end(); ++i) {sort(i->begin(), i->end(), cmp_x1);}
		if (!is_tile) {timer.end();} // use a single timer for tile mode

		if (params.flatten_mesh && !use_city_plots) { // not needed for city plots, which are already flat
			timer_t timer("Gen Building Zvals", !is_tile);
			bool const do_flatten(using_tiled_terrain_hmap_tex());

#pragma omp parallel for schedule(static,1) if (!is_tile)
			for (int i = 0; i < (int)buildings.size(); ++i) {
				building_t &b(buildings[i]);

				if (do_flatten) {
					//assert(!b.is_rotated()); // too strong?
					flatten_hmap_region(b.bcube); // flatten the mesh under the bcube to a height of mesh_zval
				}
				else { // extend building bottom downward to min mesh height
					float &zmin(b.bcube.d[2][0]); // Note: grid bcube z0 value won't be correct, but will be fixed conservatively below
					float const zmin0(zmin);
					unsigned num_below(0);
					
					for (int d = 0; d < 4; ++d) {
						float const zval(get_exact_zval(b.bcube.d[0][d&1]+xlate.x, b.bcube.d[1][d>>1]+xlate.y)); // approximate for rotated buildings
						min_eq(zmin, zval);
						num_below += (zval < def_water_level);
					}
					max_eq(zmin, def_water_level); // don't go below the water
					float const max_dz(b.get_material().max_delta_z);
					if (num_below > 2 || // more than 2 corners underwater
						(max_dz > 0.0 && (zmin0 - zmin) > max_dz)) // too steep of a slope
					{
						b.bcube.set_to_zeros();
						++num_skip;
					}
					else if (!b.parts.empty()) {
						b.parts.back().d[2][0] = b.bcube.d[2][0]; // update base z1
						assert(b.parts.back().dz() > 0.0);
					}
				}
			} // for i
			if (do_flatten) { // use conservative zmin for grid
				for (auto i = grid.begin(); i != grid.end(); ++i) {i->bcube.d[2][0] = def_water_level;}
			}
		} // if flatten_mesh
		{ // open a scope
			timer_t timer2("Gen Building Geometry", !is_tile);
#pragma omp parallel for schedule(static,1) if (!is_tile)
			for (int i = 0; i < (int)buildings.size(); ++i) {buildings[i].gen_geometry(i, 1337*i+rseed);}
		} // close the scope
		if (!is_tile) {
			cout << "WM: " << world_mode << " MCF: " << max_consec_fail << " Buildings: " << params.num_place << " / " << num_tries << " / " << num_gen
				 << " / " << buildings.size() << " / " << (buildings.size() - num_skip) << endl;
		}
		build_grid_by_tile(is_tile);
		create_vbos(is_tile);
	}

	static void multi_draw(bool shadow_only, vector3d const &xlate, vector<building_creator_t *> const &bcs) {
		if (bcs.empty()) return;
		//timer_t timer(string("Draw Buildings") + (shadow_only ? " Shadow" : "")); // 0.57ms (2.6ms with glFinish())
		bool have_windows(0), have_wind_lights(0);
		unsigned max_draw_ix(0);
		bool const night(is_night(WIND_LIGHT_ON_RAND));

		for (auto i = bcs.begin(); i != bcs.end(); ++i) {
			assert(*i);
			if (night) {(*i)->ensure_window_lights_vbos();}
			have_windows     |= !(*i)->building_draw_windows.empty();
			have_wind_lights |= !(*i)->building_draw_wind_lights.empty();
			max_eq(max_draw_ix,  (*i)->building_draw_vbo.get_num_draw_blocks());
		}
		point const camera(get_camera_pos()), camera_xlated(camera - xlate);
		int const use_bmap(global_building_params.has_normal_map);
		bool const use_tt_smap(check_tile_smap(shadow_only) && (light_valid_and_enabled(0) || light_valid_and_enabled(1))); // check for sun or moon
		shader_t s;
		fgPushMatrix();
		translate_to(xlate);
		glDepthFunc(GL_LEQUAL);

		// main/batched draw pass
		if (shadow_only) {s.begin_color_only_shader();} // really don't even need colors
		else {
			bool const v(world_mode == WMODE_GROUND), indir(v), dlights(v), use_smap(v);
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, dlights, 0, 0, (use_smap ? 2 : 1), use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw.init_draw_frame();}
		}
		for (unsigned ix = 0; ix < max_draw_ix; ++ix) {
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_vbo.draw_block(ix, shadow_only);}
		}
		if (!shadow_only && have_windows) { // draw windows
			enable_blend();
			glDepthMask(GL_FALSE); // disable depth writing
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_windows.draw(0);} // draw windows on top of other buildings
			glDepthMask(GL_TRUE); // re-enable depth writing
			disable_blend();
		}
		s.end_shader();

		// post-pass to render buildings in nearby tiles that have shadow maps
		if (use_tt_smap) {
			//timer_t timer2("Draw Buildings Smap"); // 0.3
			city_shader_setup(s, 1, 1, use_bmap, 0.0); // use_smap=1, use_dlights=1, min_alpha=0.0 to avoid alpha test
			float const draw_dist(get_tile_smap_dist() + 0.5*(X_SCENE_SIZE + Y_SCENE_SIZE));
			glDepthMask(GL_FALSE); // disable depth writing

			for (auto i = bcs.begin(); i != bcs.end(); ++i) {
				for (auto g = (*i)->grid_by_tile.begin(); g != (*i)->grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
					if (!g->bcube.closest_dist_less_than(camera_xlated, draw_dist)) continue; // too far
					point const pos(g->bcube.get_cube_center() + xlate);
					if (!camera_pdu.sphere_and_cube_visible_test(pos, g->bcube.get_bsphere_radius(), (g->bcube + xlate))) continue; // VFC
					if (!try_bind_tile_smap_at_point(pos, s)) continue; // no shadow maps - not drawn in this pass
					unsigned const tile_id(g - (*i)->grid_by_tile.begin());
					(*i)->building_draw_vbo.draw_tile(tile_id);

					if (!(*i)->building_draw_windows.empty()) {
						enable_blend();
						(*i)->building_draw_windows.draw_tile(tile_id); // draw windows on top of other buildings
						disable_blend();
					}
				} // for g
			}
			s.end_shader();
			glDepthMask(GL_TRUE); // re-enable depth writing
		}
		if (!shadow_only && night && have_wind_lights) { // add night time random lights in windows
			enable_blend();
			glDepthMask(GL_FALSE); // disable depth writing
			float const low_v(0.5 - WIND_LIGHT_ON_RAND), high_v(0.5 + WIND_LIGHT_ON_RAND), lit_thresh_mult(1.0 + 2.0*CLIP_TO_01((light_factor - low_v)/(high_v - low_v)));
			s.end_shader();
			s.set_vert_shader("window_lights");
			s.set_frag_shader("linear_fog.part+window_lights");
			s.set_prefix("#define FOG_FADE_TO_TRANSPARENT", 1);
			setup_tt_fog_pre(s);
			s.begin_shader();
			s.add_uniform_float("lit_thresh_mult", lit_thresh_mult); // gradual transition of lit window probability around sunset
			setup_tt_fog_post(s);
			for (auto i = bcs.begin(); i != bcs.end(); ++i) {(*i)->building_draw_wind_lights.draw(0);} // add bloom?
			glDepthMask(GL_TRUE); // re-enable depth writing
			disable_blend();
		}
		glDepthFunc(GL_LESS);
		fgPopMatrix();
	}

	void get_all_window_verts(building_draw_t &bdraw, bool light_pass) {
		bdraw.clear();

		for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
			bdraw.cur_tile_id = (g - grid_by_tile.begin());
			for (auto i = g->ixs.begin(); i != g->ixs.end(); ++i) {get_building(*i).get_all_drawn_window_verts(bdraw, light_pass);}
		}
		bdraw.finalize(grid_by_tile.size());
	}
	void get_all_drawn_verts() { // Note: non-const; building_draw is modified
#pragma omp parallel for schedule(static) num_threads(2)
		for (int pass = 0; pass < 2; ++pass) { // parallel loop doesn't help much because pass 0 takes most of the time
			if (pass == 0) { // main pass
				building_draw_vbo.clear();

				for (auto g = grid_by_tile.begin(); g != grid_by_tile.end(); ++g) { // Note: all grids should be nonempty
					building_draw_vbo.cur_tile_id = (g - grid_by_tile.begin());
					for (auto i = g->ixs.begin(); i != g->ixs.end(); ++i) {get_building(*i).get_all_drawn_verts(building_draw_vbo);}
				}
				building_draw_vbo.finalize(grid_by_tile.size());
			}
			else if (pass == 1) { // windows pass
				get_all_window_verts(building_draw_windows, 0);
				if (is_night(WIND_LIGHT_ON_RAND)) {get_all_window_verts(building_draw_wind_lights, 1);} // only generate window verts at night
			}
		} // for pass
	}
	void create_vbos(bool is_tile) { // Note: non-const; building_draw is modified
		building_window_gen.check_windows_texture();
		timer_t timer("Create Building VBOs", !is_tile);
		get_all_drawn_verts();
		unsigned const num_verts(building_draw_vbo.num_verts()), num_tris(building_draw_vbo.num_tris());
		if (!is_tile) {cout << "Building verts: " << num_verts << ", tris: " << num_tris << ", mem: " << num_verts*sizeof(vert_norm_comp_tc_color) << endl;}
		building_draw_vbo.upload_to_vbos();
		building_draw_windows.upload_to_vbos();
		building_draw_wind_lights.upload_to_vbos(); // Note: may be empty if not night time
	}
	void ensure_window_lights_vbos() {
		if (!building_draw_wind_lights.empty()) return; // already calculated
		building_window_gen.check_windows_texture();
		get_all_window_verts(building_draw_wind_lights, 1);
		building_draw_wind_lights.upload_to_vbos();
	}
	void clear_vbos() {
		building_draw.clear_vbos();
		building_draw_vbo.clear_vbos();
		building_draw_windows.clear_vbos();
		building_draw_wind_lights.clear_vbos();
	}

	bool check_sphere_coll(point &pos, point const &p_last, float radius, bool xy_only=0, vector3d *cnorm=nullptr) const {
		if (empty()) return 0;
		vector3d const xlate(get_camera_coord_space_xlate());
		cube_t bcube; bcube.set_from_sphere((pos - xlate), radius);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		float const dist(p2p_dist(pos, p_last));
		vector<point> points; // reused across calls

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.ixs.empty()) continue; // skip empty grid
				if (!sphere_cube_intersect(pos, (radius + dist), (ge.bcube + xlate))) continue; // Note: makes little difference

				// Note: assumes buildings are separated so that only one sphere collision can occur
				for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) {
					building_t const &building(get_building(*b));
					if (!building.bcube.intersects_xy(bcube)) continue;
					if (building.check_sphere_coll(pos, p_last, xlate, radius, xy_only, points, cnorm)) return 1;
				}
			} // for x
		} // for y
		return 0;
	}

	unsigned check_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix, bool ret_any_pt) const {
		if (empty()) return 0;
		bool const is_vertical(p1.x == p2.x && p1.y == p2.y); // vertical lines can only intersect one building
		vector3d const xlate(get_camera_coord_space_xlate());
		point const p1x(p1 - xlate);
		cube_t bcube(p1x, p2-xlate);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		point end_pos(p2);
		unsigned coll(0); // 0=none, 1=side, 2=roof
		vector<point> points; // reused across calls
		t = 1.0; // start at end point

		// for now, just do a slow iteration over every grid element within the line's bbox in XY
		// Note: should probably iterate over the grid in XY order from the start to the end of the line, or better yet use a line drawing algorithm
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (ge.ixs.empty()) continue; // skip empty grid
				if (!check_line_clip(p1x, (end_pos - xlate), ge.bcube.d)) continue; // no intersection - skip this grid

				for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) { // Note: okay to check the same building more than once
					building_t const &building(get_building(*b));
					if (!building.bcube.intersects(bcube)) continue;
					float t_new(t);
					unsigned const ret(building.check_line_coll(p1, p2, xlate, t_new, points, 0, ret_any_pt));

					if (ret && t_new <= t) { // closer hit pos, update state
						t       = t_new;
						hit_bix = *b;
						coll    = ret;
						end_pos = p1 + t*(p2 - p1);
						if (ret_any_pt || is_vertical) return coll;
					}
				}
			} // for x
		} // for y
		return coll;
	}

	// Note: we can get building_id by calling check_ped_coll() or get_building_bcube_at_pos()
	bool check_line_coll_building(point const &p1, point const &p2, unsigned building_id) const { // Note: not thread safe sue to static points
		assert(building_id < buildings.size());
		static vector<point> points; // reused across calls
		float t_new(1.0);
		return buildings[building_id].check_line_coll(p1, p2, zero_vector, t_new, points, 0, 1);
	}

	int get_building_bcube_contains_pos(point const &pos) { // Note: not thread safe sue to static points
		if (empty()) return -1;
		unsigned const gix(get_grid_ix(pos));
		grid_elem_t const &ge(grid[gix]);
		if (ge.ixs.empty() || !ge.bcube.contains_pt(pos)) return -1; // skip empty or non-containing grid
		static vector<point> points; // reused across calls

		for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) {
			if (get_building(*b).bcube.contains_pt(pos)) {return *b;} // found
		}
		return -1;
	}

	bool check_ped_coll(point const &pos, float radius, unsigned plot_id, unsigned &building_id) const { // Note: not thread safe sue to static points
		if (empty()) return 0;
		assert(plot_id < bix_by_plot.size());
		vector<unsigned> const &bixes(bix_by_plot[plot_id]); // should be populated in gen()
		if (bixes.empty()) return 0;
		cube_t bcube; bcube.set_from_sphere(pos, radius);
		static vector<point> points; // reused across calls

		// Note: assumes buildings are separated so that only one ped collision can occur
		for (auto b = bixes.begin(); b != bixes.end(); ++b) {
			building_t const &building(get_building(*b));
			if (building.bcube.x1() > bcube.x2()) break; // no further buildings can intersect (sorted by x1)
			if (!building.bcube.intersects_xy(bcube)) continue;
			if (building.check_point_or_cylin_contained(pos, 2.0*radius, points)) {building_id = *b; return 1;} // double the radius value to add padding to account for inaccuracy
		}
		return 0;
	}
	bool select_building_in_plot(unsigned plot_id, unsigned rand_val, unsigned &building_id) const {
		if (bix_by_plot.empty()) return 0; // not setup / no buildings
		assert(plot_id < bix_by_plot.size());
		vector<unsigned> const &bixes(bix_by_plot[plot_id]);
		if (bixes.empty()) return 0;
		building_id = bixes[rand_val % bixes.size()];
		return 1;
	}

	void get_overlapping_bcubes(cube_t const &xy_range, vect_cube_t &bcubes) const { // Note: called on init, don't need to use get_camera_coord_space_xlate()
		if (empty()) return; // nothing to do
		unsigned ixr[2][2];
		get_grid_range(xy_range, ixr);

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (!xy_range.intersects_xy(ge.bcube)) continue;

				for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) {
					building_t const &building(get_building(*b));
					if (!xy_range.intersects_xy(building.bcube)) continue;
					cube_t shared(xy_range);
					shared.intersect_with_cube(building.bcube);
					if (get_grid_ix(shared.get_llc()) == y*grid_sz + x) {bcubes.push_back(building.bcube);} // add only if in home grid (to avoid duplicates)
				}
			} // for x
		} // for y
	}

	void get_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state) const {
		state.init(pdu.pos, get_camera_coord_space_xlate());
		
		for (auto g = grid.begin(); g != grid.end(); ++g) {
			if (g->ixs.empty()) continue;
			point const pos(g->bcube.get_cube_center() + state.xlate);
			if (!pdu.sphere_and_cube_visible_test(pos, g->bcube.get_bsphere_radius(), (g->bcube + state.xlate))) continue; // VFC
			
			for (auto i = g->ixs.begin(); i != g->ixs.end(); ++i) {
				if (pdu.cube_visible(buildings[*i].bcube + state.xlate)) {state.building_ids.push_back(*i);}
			}
		}
	}
	bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state) const {
		for (vector<unsigned>::const_iterator b = state.building_ids.begin(); b != state.building_ids.end(); ++b) {
			building_t const &building(get_building(*b));
			bool occluded(1);

			for (unsigned i = 0; i < npts; ++i) {
				float t(1.0); // start at end of line
				if (!building.check_line_coll(state.pos, pts[i], state.xlate, t, state.temp_points, 1)) {occluded = 0; break;}
			}
			if (occluded) return 1;
		} // for b
		return 0;
	}
}; // building_creator_t


class building_tiles_t {
	map<pair<int, int>, building_creator_t> tiles; // key is {x, y} pair
public:
	bool empty() const {return tiles.empty();}

	bool create_tile(int x, int y) {
		auto it(tiles.find(make_pair(x, y)));
		if (it != tiles.end()) return 0; // already exists
		//cout << "Create building tile " << x << "," << y << ", tiles: " << tiles.size() << endl; // 299 tiles
		building_creator_t &bc(tiles[make_pair(x, y)]); // insert it
		assert(bc.empty());
		cube_t bcube(all_zeros);
		bcube.x1() = get_xval(x*MESH_X_SIZE);
		bcube.y1() = get_yval(y*MESH_Y_SIZE);
		bcube.x2() = get_xval((x+1)*MESH_X_SIZE);
		bcube.y2() = get_yval((y+1)*MESH_Y_SIZE);
		global_building_params.set_pos_range(bcube);
		// TODO: Optimize Gen + Optimize Draw
		int const rseed(x + (y << 16) + 12345); // should not be zero
		bc.gen(global_building_params, 0, 0, 1, rseed);
		global_building_params.restore_prev_pos_range();
		return 1;
	}
	bool remove_tile(int x, int y) {
		auto it(tiles.find(make_pair(x, y)));
		//cout << "Remove building tile " << x << "," << y << ", tiles: " << tiles.size() << endl;
		if (it == tiles.end()) return 0; // not found
		it->second.clear(); // free VBOs/VAOs
		tiles.erase(it);
		return 1;
	}
	void clear_vbos() {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {i->second.clear_vbos();}
	}
	void clear() {
		clear_vbos();
		tiles.clear();
	}
	vector3d get_max_extent() const {
		vector3d max_extent(zero_vector);
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {max_extent = max_extent.max(i->second.get_max_extent());}
		return max_extent;
	}
	bool check_sphere_coll(point &pos, point const &p_last, float radius, bool xy_only=0, vector3d *cnorm=nullptr) const {
		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			if (i->second.check_sphere_coll(pos, p_last, radius, xy_only, cnorm)) return 1;
		}
		return 0;
	}
	void add_drawn(vector3d const &xlate, vector<building_creator_t *> &bcs) {
		float const draw_dist(get_draw_tile_dist());
		point const camera(get_camera_pos() - xlate);

		for (auto i = tiles.begin(); i != tiles.end(); ++i) {
			//if (!i->second.get_bcube().closest_dist_xy_less_than(camera, draw_dist)) continue; // distance test (conservative)
			if (!dist_xy_less_than(camera, i->second.get_bcube().get_cube_center(), draw_dist)) continue; // distance test (aggressive)
			if (i->second.is_visible(xlate)) {bcs.push_back(&i->second);}
		}
	}
}; // end building_tiles_t


building_creator_t building_creator, building_creator_city;
building_tiles_t building_tiles;

void create_buildings_tile(int x, int y) {
	if (global_building_params.gen_inf_buildings()) {building_tiles.create_tile(x, y);}
}
void remove_buildings_tile(int x, int y) {
	if (global_building_params.gen_inf_buildings()) {building_tiles.remove_tile(x, y);}
}

vector3d get_tt_xlate_val() {return ((world_mode == WMODE_INF_TERRAIN) ? vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0) : zero_vector);}

void gen_buildings() {
	update_sun_and_moon(); // need to update light_factor from sun to know if we need to generate window light geometry

	if (have_cities()) {
		building_creator_city.gen(global_building_params, 1, 0, 0); // city buildings
		global_building_params.restore_prev_pos_range(); // hack to undo clip to city bounds to allow buildings to extend further out
		building_creator.gen     (global_building_params, 0, 1, 0); // non-city secondary buildings
	} else {building_creator.gen (global_building_params, 0, 0, 0);} // mixed buildings
}
void draw_buildings(bool shadow_only, vector3d const &xlate) {
	if (world_mode != WMODE_INF_TERRAIN) {building_tiles.clear();}
	vector<building_creator_t *> bcs;
	if (building_creator_city.is_visible(xlate)) {bcs.push_back(&building_creator_city);}
	if (building_creator     .is_visible(xlate)) {bcs.push_back(&building_creator     );}
	building_tiles.add_drawn(xlate, bcs);
	building_creator_t::multi_draw(shadow_only, xlate, bcs);
}
bool proc_buildings_sphere_coll(point &pos, point const &p_int, float radius, bool xy_only, vector3d *cnorm) {
	return (building_creator_city.check_sphere_coll(pos, p_int, radius, xy_only, cnorm) ||
		         building_creator.check_sphere_coll(pos, p_int, radius, xy_only, cnorm) ||
		           building_tiles.check_sphere_coll(pos, p_int, radius, xy_only, cnorm));
}
bool check_buildings_sphere_coll(point const &pos, float radius, bool apply_tt_xlate, bool xy_only) {
	point center(pos);
	if (apply_tt_xlate) {center += get_tt_xlate_val();} // apply xlate for all static objects - not the camera
	return proc_buildings_sphere_coll(center, pos, radius, xy_only, nullptr);
}
bool check_buildings_point_coll(point const &pos, bool apply_tt_xlate, bool xy_only) {
	return check_buildings_sphere_coll(pos, 0.0, apply_tt_xlate, xy_only);
}
unsigned check_buildings_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix, bool apply_tt_xlate, bool ret_any_pt) { // for line_intersect_city()
	vector3d const xlate(apply_tt_xlate ? get_tt_xlate_val() : zero_vector);
	return (building_creator_city.check_line_coll(p1+xlate, p2+xlate, t, hit_bix, ret_any_pt) ||
		         building_creator.check_line_coll(p1+xlate, p2+xlate, t, hit_bix, ret_any_pt));
}
bool get_buildings_line_hit_color(point const &p1, point const &p2, colorRGBA &color) {
	// only check city buildings in TT mode (too slow to check others)
	return ((world_mode == WMODE_GROUND) ? building_creator : building_creator_city).get_building_hit_color(p1, p2, color);
}
bool have_buildings() {return (!building_creator.empty() || !building_creator_city.empty() || !building_tiles.empty());} // for postproce effects

vector3d get_buildings_max_extent() { // used for TT shadow bounds + map mode
	return building_creator.get_max_extent().max(building_creator_city.get_max_extent()).max(building_tiles.get_max_extent());
}
void clear_building_vbos() {
	building_creator.clear_vbos();
	building_creator_city.clear_vbos();
	building_tiles.clear_vbos();
}

// city interface
void set_buildings_pos_range(cube_t const &pos_range) {global_building_params.set_pos_range(pos_range);}
void get_building_bcubes(cube_t const &xy_range, vect_cube_t &bcubes) {building_creator_city.get_overlapping_bcubes(xy_range, bcubes);} // Note: no xlate applied
// cars + peds
void get_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state) {building_creator_city.get_occluders(pdu, state);}
bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state) {return building_creator_city.check_pts_occluded(pts, npts, state);}
// used for pedestrians
cube_t get_building_bcube(unsigned building_id) {return building_creator_city.get_building_bcube(building_id);}
bool check_line_coll_building(point const &p1, point const &p2, unsigned building_id) {return building_creator_city.check_line_coll_building(p1, p2, building_id);}
int get_building_bcube_contains_pos(point const &pos) {return building_creator_city.get_building_bcube_contains_pos(pos);}
bool check_buildings_ped_coll(point const &pos, float radius, unsigned plot_id, unsigned &building_id) {return building_creator_city.check_ped_coll(pos, radius, plot_id, building_id);}
bool select_building_in_plot(unsigned plot_id, unsigned rand_val, unsigned &building_id) {return building_creator_city.select_building_in_plot(plot_id, rand_val, building_id);}

