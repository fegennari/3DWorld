// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "file_utils.h"

using std::string;

unsigned const MAX_CYLIN_SIDES = 36;

extern int rand_gen_index, display_mode;
extern float shadow_map_pcf_offset, cobj_z_bias;

// TODO:
// Multilevel cylinders and N-gons shapes?
// Texture alignment for windows

struct tid_nm_pair_t {

	int tid, nm_tid; // Note: assumes each tid has only one nm_tid
	float tscale_x, tscale_y;

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale_x(1.0), tscale_y(1.0) {}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}
	bool operator==(tid_nm_pair_t const &t) const {return (tid == t.tid && nm_tid == t.nm_tid && tscale_x == t.tscale_x && tscale_y == t.tscale_y);}
	colorRGBA get_avg_color() const {return texture_color(tid);}

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

	unsigned min_levels, max_levels, min_sides, max_sides;
	float place_radius, max_delta_z, max_rot_angle, min_level_height, min_alt, max_alt;
	float split_prob, cube_prob, round_prob, asf_prob, min_fsa, max_fsa, min_asf, max_asf;
	cube_t pos_range, sz_range; // pos_range z is unused?
	color_range_t side_color, roof_color;

	building_mat_t() : min_levels(1), max_levels(1), min_sides(4), max_sides(4), place_radius(0.0), max_delta_z(0.0), max_rot_angle(0.0),
		min_level_height(0.0), min_alt(-1000), max_alt(1000), split_prob(0.0), cube_prob(1.0), round_prob(0.0), asf_prob(0.0),
		min_fsa(0.0), max_fsa(0.0), min_asf(0.0), max_asf(0.0), pos_range(-100,100,-100,100,0,0), sz_range(1,1,1,1,1,1) {}
	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0);}

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
};

struct building_params_t {

	bool flatten_mesh, has_normal_map, tex_mirror, tex_inv_y;
	unsigned num_place, num_tries, cur_prob;
	float ao_factor;
	vector3d range_translate; // used as a temporary to add to material pos_range
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix;

	building_params_t(unsigned num=0) : flatten_mesh(0), has_normal_map(0), tex_mirror(0), tex_inv_y(0),
		num_place(num), num_tries(10), cur_prob(1), ao_factor(0.0), range_translate(zero_vector) {}
	int get_wrap_mir() const {return (tex_mirror ? 2 : 1);}
	void add_cur_mat() {
		unsigned const mat_ix(materials.size());
		for (unsigned n = 0; n < cur_prob; ++n) {mat_gen_ix.push_back(mat_ix);} // add more references to this mat for higher probability
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
	unsigned choose_rand_mat(rand_gen_t &rgen) const {
		assert(!mat_gen_ix.empty());
		return mat_gen_ix[rgen.rand()%mat_gen_ix.size()];
	}
};

building_params_t global_building_params;

void buildings_file_err(string const &str, int &error) {
	cout << "Error reading buildings config option " << str << "." << endl;
	error = 1;
}

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
	else if (str == "side_nm_tid") {
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
	bool half_offset;
	float rot_sin, rot_cos, flat_side_amt, alt_step_factor; // rotation in XY plane, around Z (up) axis
	//float roof_recess;

	building_geom_t() : num_sides(4), half_offset(0), rot_sin(0.0), rot_cos(1.0), flat_side_amt(0.0), alt_step_factor(0.0) {}
};

struct building_t : public building_geom_t {

	unsigned mat_ix;
	colorRGBA side_color, roof_color;
	cube_t bcube;
	vector<cube_t> parts;
	mutable unsigned cur_draw_ix;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_), side_color(WHITE), roof_color(WHITE), cur_draw_ix(0) {bcube.set_to_zeros();}
	bool is_valid  () const {return !bcube.is_all_zeros();}
	bool is_rotated() const {return (rot_sin != 0.0);}
	bool is_cube()    const {return (num_sides == 4);}
	bool use_cylinder_coll() const {return (num_sides > 8 && flat_side_amt == 0.0);} // use cylinder collision if not a cube, triangle, octagon, etc. (approximate)
	colorRGBA get_avg_side_color() const {return side_color.modulate_with(get_material().side_tex.get_avg_color());}
	colorRGBA get_avg_roof_color() const {return roof_color.modulate_with(get_material().roof_tex.get_avg_color());}
	building_mat_t const &get_material() const {return global_building_params.get_material(mat_ix);}
	void gen_rotation(rand_gen_t &rgen);
	bool check_bcube_overlap_xy(building_t const &b, float expand) const {
		return (check_bcube_overlap_xy_one_dir(b, expand) || b.check_bcube_overlap_xy_one_dir(*this, expand));
	}
	bool check_sphere_coll(point const &pos, float radius, bool xy_only=0) const {
		point pos2(pos);
		return check_sphere_coll(pos2, pos, zero_vector, radius, xy_only);
	}
	bool check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only=0) const;
	unsigned check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t) const;
	void gen_geometry(unsigned ix);
	void draw(shader_t &s, bool shadow_only, float far_clip, vector3d const &xlate, building_draw_t &bdraw, unsigned draw_ix) const;
private:
	bool check_bcube_overlap_xy_one_dir(building_t const &b, float expand) const;
	void split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen);
};


void do_xy_rotate(float rot_sin, float rot_cos, point const &center, point &pos) {
	float const x(pos.x - center.x), y(pos.y - center.y); // translate to center
	pos.x = x*rot_cos - y*rot_sin + center.x;
	pos.y = y*rot_cos + x*rot_sin + center.y;
}


#define EMIT_VERTEX() \
	vert.v = pt*sz + llc; \
	vert.t[ st] = tscale[ st]*vert.v[d]; \
	vert.t[!st] = tscale[!st]*vert.v[i]; \
	if (apply_ao) {vert.copy_color(cw[pt.z == 1]);} \
	if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, center, vert.v);} \
	verts.push_back(vert);

class building_draw_t {

	struct draw_block_t {
		tid_nm_pair_t tex;
		vector<vert_norm_comp_tc_color> quad_verts, tri_verts;

		void draw_geom(bool shadow_only, int force_tid=-1) {
			if (empty()) return;
			if (force_tid >= 0) {select_texture(force_tid); select_multitex(FLAT_NMAP_TEX, 5);} // no normal map
			else if (!shadow_only) {tex.set_gl();}
			draw_quad_verts_as_tris(quad_verts);
			draw_verts(tri_verts, GL_TRIANGLES);
		}
		void draw_and_clear(bool shadow_only, int force_tid=-1) {draw_geom(shadow_only, force_tid); clear();}
		void clear() {quad_verts.clear(); tri_verts.clear();}
		bool empty() const {return (quad_verts.empty() && tri_verts.empty());}
	};
	vector<draw_block_t> to_draw, pend_draw; // one per texture, assumes tids are dense

	vector<vert_norm_comp_tc_color> &get_verts(tid_nm_pair_t const &tex, bool quads_or_tris=0) { // default is quads
		unsigned const ix((tex.tid >= 0) ? (tex.tid+1) : 0);
		if (ix >= to_draw.size()) {to_draw.resize(ix+1);}
		draw_block_t &block(to_draw[ix]);
		if (block.empty()) {block.tex = tex;} // copy material first time
		else {assert(block.tex.nm_tid == tex.nm_tid);} // else normal maps must agree
		return (quads_or_tris ? block.tri_verts : block.quad_verts);
	}
	static void setup_ao_color(colorRGBA const &color, cube_t const &bcube, float z1, float z2, color_wrapper cw[2], vert_norm_comp_tc_color &vert) {
		if (global_building_params.ao_factor > 0.0) {
			float const dz_mult(global_building_params.ao_factor/bcube.get_dz());
			UNROLL_2X(cw[i_].set_c4(color*((1.0 - global_building_params.ao_factor) + dz_mult*((i_ ? z2 : z1) - bcube.d[2][0])));)
		} else {vert.set_c4(color);} // color is shared across all verts
	}
	vector<vector3d> normals; // reused across add_cylinder() calls
	vector<point> points; // reused across calc_poly_pts() calls
	point cur_camera_pos;

public:
	building_draw_t() : cur_camera_pos(zero_vector) {}
	void init_draw_frame() {cur_camera_pos = get_camera_pos();} // capture camera pos during non-shadow pass to use for shadow pass

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
		float sin_s(0.0), cos_s(1.0); // start at 0 - more efficient
		if (bg.half_offset) {sin_s = sin(0.5*css); cos_s = cos(0.5*css);} // for cube
		nv.resize(ndiv);

		for (unsigned S = 0; S < ndiv; ++S) { // build normals table
			bool const d(S&1);
			float const s(sin_s), c(cos_s);
			nv[S].assign(s, c, 0.0);
			sin_s = s*cos_ds[d] + c*sin_ds[d];
			cos_s = c*cos_ds[d] - s*sin_ds[d];
		}
	}
	static void calc_poly_pts(building_geom_t const &bg, cube_t const &bcube, vector<point> &pts, unsigned ndiv, float expand=0.0) {
		calc_normals(bg, pts, ndiv);
		vector3d const sz(bcube.get_size());
		point const cc(bcube.get_cube_center());
		float const rscale(0.5), rx(rscale*sz.x + expand), ry(rscale*sz.y + expand); // expand polygon by sphere radius
		for (unsigned i = 0; i < ndiv; ++i) {pts[i].assign((cc.x + rx*pts[i].x), (cc.y + ry*pts[i].y), 0.0);} // convert normals to points
	}

	void add_cylinder(building_geom_t const &bg, point const &pos, point const &rot_center, float height, float rx, float ry, point const &xlate,
		cube_t const &bcube, tid_nm_pair_t const &tex, colorRGBA const &color, bool shadow_only, vector3d const &view_dir, unsigned dim_mask)
	{
		unsigned ndiv(bg.num_sides);
		assert(ndiv >= 3);
		bool const smooth_normals(ndiv >= 16); // cylinder vs. N-gon
		
		if (ndiv > 4 && bg.flat_side_amt == 0.0 && bg.alt_step_factor == 0.0) {
			float const dist(max(p2p_dist(cur_camera_pos, (pos + xlate)), 0.001f));
			ndiv = max(min(ndiv, unsigned(1000.0*max(rx, ry)/dist)), 3U); // LOD if not flat sides: use at least 3 sides
		}
		float const ndiv_inv(1.0/ndiv), z_top(pos.z + height), tscale_x(2.0*tex.tscale_x), tscale_y(2.0*tex.tscale_y); // adjust for local vs. global space change
		bool const apply_ao(!shadow_only && global_building_params.ao_factor > 0.0);
		vert_norm_comp_tc_color vert;
		color_wrapper cw[2];
		setup_ao_color(color, bcube, pos.z, z_top, cw, vert);
		float tex_pos[2] = {0.0, 1.0};
		calc_normals(bg, normals, ndiv);

		if (!shadow_only) {
			float const dz_inv(1.0/bcube.get_dz());
			UNROLL_2X(tex_pos[i_] = dz_inv*((i_ ? z_top : pos.z) - bcube.d[2][0]);)
		}
		if (dim_mask & 3) { // draw sides
			auto &verts(get_verts(tex)); // Note: cubes are drawn with quads, so we want to emit quads here

			if (shadow_only) {
				for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data quads
					for (unsigned d = 0; d < 2; ++d) {
						vector3d const &n(normals[d ? ((S+1)%ndiv) : S]);
						vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), 0.0);
						if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, vert.v);}

						for (unsigned e = 0; e < 2; ++e) {
							vert.v.z = ((d^e) ? z_top : pos.z);
							verts.push_back(vert);
						}
					} // for d
				} // for S
			}
			else {
				float tot_perim(0.0), tot_perim_inv(0.0), cur_perim[2] = {0.0, 0.0};
				for (unsigned S = 0; S < ndiv; ++S) {tot_perim += p2p_dist(normals[S], normals[(S+1)%ndiv]);}
				tot_perim_inv = 1.0/tot_perim;
				
				for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data quads
					vector3d const &n1(normals[S]), &n2(normals[(S+1)%ndiv]);
					cur_perim[0]  = cur_perim[1];
					cur_perim[1] += p2p_dist(n1, n2);
					vector3d normal(n1 + n2); normal.x *= ry; normal.y *= rx; // average the two vertex normals for the flat face normal
					if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, all_zeros, normal);}
					if (view_dir.x*normal.x + view_dir.y*normal.y > 0.0) continue; // back facing
					if (!smooth_normals) {vert.set_norm(normal.get_norm());}

					for (unsigned d = 0; d < 2; ++d) {
						vector3d const &n(d ? n2 : n1);
						vert.t[0] = tscale_x*cur_perim[d]*tot_perim_inv; // texture_scale should be a multiple of 1.0
					
						if (smooth_normals) {
							vector3d normal(n); normal.x *= ry; normal.y *= rx; // scale normal by radius (swapped)
							if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, all_zeros, normal);}
							vert.set_norm(normal.get_norm());
						}
						vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), 0.0);
						if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, vert.v);}

						for (unsigned e = 0; e < 2; ++e) {
							vert.v.z = ((d^e) ? z_top : pos.z);
							vert.t[1] = tscale_y*tex_pos[d^e];
							if (apply_ao) {vert.copy_color(cw[d^e]);}
							verts.push_back(vert);
						}
					} // for d
				} // for S
			}
		} // end draw sides
		if (dim_mask & 4) { // draw end(s) / roof
			auto &tri_verts(get_verts(tex, 1));
			
			for (unsigned d = 0; d < 2; ++d) { // bottom, top
				if ((view_dir.z < 0.0) ^ d) continue; // back facing
				float const zval(d ? z_top : pos.z);
				vert.set_norm(d ? plus_z : -plus_z);
				if (apply_ao) {vert.copy_color(cw[d]);}
				vert_norm_comp_tc_color center(vert);
				center.t[0] = center.t[1] = 0.0; // center of texture space for this disk
				center.v = pos;
				if (d) {center.v.z += height;}
				if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, center.v);}

				for (unsigned S = 0; S < ndiv; ++S) { // generate vertex data triangles
					tri_verts.push_back(center);

					for (unsigned e = 0; e < 2; ++e) {
						if (S > 0 && e == 0) {tri_verts.push_back(tri_verts[tri_verts.size()-2]); continue;} // reuse prev vertex
						vector3d const &n(normals[(S+e)%ndiv]);
						vert.v.assign((pos.x + rx*n.x), (pos.y + ry*n.y), center.v.z);
						if (!shadow_only) {vert.t[0] = tscale_x*n[0]; vert.t[1] = tscale_y*n[1];}
						if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, rot_center, vert.v);}
						tri_verts.push_back(vert);
					}
				} // for S
			} // for d
		} // end draw end(s)
	}

	void add_section(building_geom_t const &bg, cube_t const &cube, point const &xlate, cube_t const &bcube,
		tid_nm_pair_t const &tex, colorRGBA const &color, bool shadow_only, vector3d const &view_dir, unsigned dim_mask)
	{
		assert(bg.num_sides >= 3); // must be nonzero volume
		point const center((bg.rot_sin == 0.0) ? all_zeros : bcube.get_cube_center()); // rotate about bounding cube / building center
		vector3d const sz(cube.get_size());

		if (bg.num_sides != 4) { // not a cube, use cylinder
			vector3d const bcube_sz(bcube.get_size());
			point const ccenter(cube.get_cube_center()), pos(ccenter.x, ccenter.y, cube.d[2][0]);
			//float const rscale(0.5*((num_sides <= 8) ? SQRT2 : 1.0)); // larger for triangles/cubes/hexagons/octagons (to ensure overlap/connectivity), smaller for cylinders
			float const rscale(0.5); // use shape contained in bcube so that bcube tests are correct, since we're not creating L/T/U shapes for this case
			add_cylinder(bg, pos, center, sz.z, rscale*sz.x, rscale*sz.y, xlate, bcube, tex, color, shadow_only, view_dir, dim_mask);
			return;
		}
		// else draw as a cube (optimized flow)
		auto &verts(get_verts(tex));
		vector3d const llc(cube.get_llc()); // move origin from center to min corner
		vert_norm_comp_tc_color vert;

		if (shadow_only) {
			for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
				unsigned const n((i+2)%3);
				if (!(dim_mask & (1<<n))) continue;
				unsigned const d((i+1)%3);
				for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
					point pt; pt[n] = j;
					//if (bg.roof_recess > 0.0 && n == 2 && j == 1) {pt.z -= bg.roof_recess*cube.get_dz();}
					for (unsigned s1 = 0; s1 < 2; ++s1) {
						pt[d] = s1;
						for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
							pt[i]  = k^j^s1^1; // need to orient the vertices differently for each side
							vert.v = pt*sz + llc;
							if (bg.rot_sin != 0.0) {do_xy_rotate(bg.rot_sin, bg.rot_cos, center, vert.v);}
							verts.push_back(vert);
						}
					}
				} // for j
			} // for i
			return;
		}
		float const tscale[2] = {2.0f*tex.tscale_x, 2.0f*tex.tscale_y}; // adjust for local vs. global space change
		bool const apply_ao(global_building_params.ao_factor > 0.0);
		color_wrapper cw[2];
		setup_ao_color(color, bcube, cube.d[2][0], cube.d[2][1], cw, vert);
		
		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const n((i+2)%3);
			if (!(dim_mask & (1<<n))) continue;
			unsigned const d((i+1)%3);
			bool const st(i&1);

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (n < 2 && bg.rot_sin != 0.0) { // XY only
					vector3d norm; norm.z = 0.0;
					if (n == 0) {norm.x =  bg.rot_cos; norm.y = bg.rot_sin;} // X
					else        {norm.x = -bg.rot_sin; norm.y = bg.rot_cos;} // Y
					if ((view_dir.x*norm.x + view_dir.y*norm.y < 0.0) ^ j) continue; // back facing
					vert.set_norm(j ? norm : -norm);
				}
				else {
					if ((view_dir[n] < 0.0) ^ j) continue; // back facing
					vert.n[i] = 0;
					vert.n[d] = 0;
					vert.n[n] = (j ? 127 : -128); // -1.0 or 1.0
				}
				point pt;
				pt[n] = j; // in direction or normal
				pt[d] = 0;
				pt[i] = !j; // need to orient the vertices differently for each side
				//if (bg.roof_recess > 0.0 && n == 2 && j == 1) {pt.z -= bg.roof_recess*cube.get_dz();}
				EMIT_VERTEX();
				pt[i] = j;
				EMIT_VERTEX();
				pt[d] = 1;
				EMIT_VERTEX();
				pt[i] = !j;
				EMIT_VERTEX();
			} // for j
		} // for i
	}
	void draw_and_clear(bool shadow_only) {
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->draw_and_clear(shadow_only);}
	}
	void begin_immediate_building() { // to be called before any add_section() calls
		pend_draw.swap(to_draw); // move current draw queue to pending queue
	}
	void end_immediate_building(bool shadow_only) { // to be matched with begin_building()
		// Note: in this case, there generally aren't more than one building of the same material within the same tile, so batching doesn't help
		draw_and_clear(shadow_only); // draw current building - sparse iteration?
		pend_draw.swap(to_draw); // restore draw queue
	}
};

building_draw_t building_draw;


void building_t::split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen) {

	// generate L, T, or U shape
	point const llc(seed_cube.get_llc()), sz(seed_cube.get_size());
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // {x,y}, {neg,pos}
	int const shape(rand()%7); // 0-6
	float const div(rgen.rand_uniform(0.3, 0.7)), s1(rgen.rand_uniform(0.2, 0.4)), s2(rgen.rand_uniform(0.6, 0.8)); // split pos in 0-1 range
	float const dpos(llc[dim] + div*sz[dim]), spos1(llc[!dim] + s1*sz[!dim]), spos2(llc[!dim] + s2*sz[!dim]); // split pos in cube space
	unsigned const start(parts.size()), num((shape == 6) ? 3 : 2);
	parts.resize(start+num, seed_cube);
	parts[start].d[dim][dir] = dpos; // full width part
	for (unsigned n = 1; n < num; ++n) {parts[start+n].d[dim][!dir] = dpos;} // partial width parts

	switch (shape) {
	case 0: case 1: case 2: case 3: // L
		parts[start+1].d[!dim][shape>>1] = ((shape&1) ? spos2 : spos1);
		break;
	case 4: case 5: // T
		parts[start+1].d[!dim][0] = spos1;
		parts[start+1].d[!dim][1] = spos2;
		break;
	case 6: // U
		parts[start+1].d[!dim][1] = spos1;
		parts[start+2].d[!dim][0] = spos2;
		break;
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
		if (i == 0) {bcube.set_from_point(corner);} else {bcube.union_with_pt(corner);}
	}
}

bool building_t::check_bcube_overlap_xy_one_dir(building_t const &b, float expand) const { // can be called before levels/splits are created

	if (expand == 0.0 && !bcube.intersects(b.bcube)) return 0;
	if (!is_rotated() && !b.is_rotated()) return 1; // above check is exact, top-level bcube check up to the caller
	point const center1(b.bcube.get_cube_center()), center2(bcube.get_cube_center());
	
	for (auto p1 = b.parts.begin(); p1 != b.parts.end(); ++p1) {
		point pts[5];
		pts[0] = p1->get_cube_center();
		cube_t c_exp(*p1);
		c_exp.expand_by(expand*p1->get_size());

		for (unsigned i = 0; i < 4; ++i) {
			pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0); // XY only
			do_xy_rotate(b.rot_sin, b.rot_cos, center1, pts[i+1]); // rotate into global space (pts[0] doesn't change)
		}
		for (unsigned i = 0; i < 5; ++i) {do_xy_rotate(-rot_sin, rot_cos, center2, pts[i]);} // inverse rotate into local coord space - negate the sine term
		
		for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) {
			for (unsigned i = 0; i < 5; ++i) {if (p2->contains_pt_xy(pts[i])) return 1;}
		}
	}
	return 0;
}

bool building_t::check_sphere_coll(point &pos, point const &p_last, vector3d const &xlate, float radius, bool xy_only) const {

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

		if (use_cylinder_coll()) {
			if ((pos2.z + radius < i->d[2][0] + xlate.z) || (pos2.z - radius > i->d[2][1] + xlate.z)) continue; // test z overlap
			point const cc(i->get_cube_center() + xlate);
			vector3d const csz(i->get_size());
#if 0 // unstable
			float const dx(pos2.x - cc.x), dy(pos2.y - cc.y), rx(0.5*csz.x + radius), ry(0.5*csz.y + radius);
			if (dx*dx/(rx*rx) + dy*dy/(ry*ry) > 1.0) continue; // approximate ellipse test
			float const mag(sqrt(dx*dx + dy*dy));
			pos2.x = cc.x + dx*rx/mag;
			pos2.y = cc.y + dy*ry/mag;
			had_coll = 1;
#else // conservative
			float const r_sum(radius + 0.5*max(csz.x, csz.y));

			if (dist_xy_less_than(pos2, cc, r_sum)) { // FIXME: conservative for ellipse case
				vector3d const dir(vector3d(pos2.x-cc.x, pos2.y-cc.y, 0.0).get_norm()); // xy dir
				UNROLL_2X(pos2[i_] = cc[i_] + dir[i_]*r_sum;)
				had_coll = 1;
			}
#endif
		}
		else if (num_sides != 4) { // triangle, hexagon, octagon, etc.
			vector<point> points;
			building_draw.calc_poly_pts(*this, (*i + xlate), points, num_sides); // without the expand
			float const dist(p2p_dist(p_last2, pos2));
			point quad_pts[4]; // quads
			bool updated(0);

			// FIXME: if the player is moving too quickly, the intersection with a side polygon may be missed,
			// which allows the player to travel through the building, but using a line intersection test from p_past2 to pos2 has other problems
			for (unsigned S = 0; S < num_sides; ++S) { // generate vertex data quads
				for (unsigned d = 0, ix = 0; d < 2; ++d) {
					point const &p(points[(S+d)%num_sides]);
					for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, i->d[2][d^e]);}
				}
				vector3d const normal(get_poly_norm(quad_pts));
				float const rdist(dot_product_ptv(normal, pos2, quad_pts[0]));
				if (rdist < 0.0 || rdist >= radius) continue; // too far or wrong side
				if (!sphere_poly_intersect(quad_pts, 4, pos2, normal, rdist, radius)) continue;
				pos2 += normal*min(rdist, dist); // calculate intersection point, adjust outward by min of distance and step size (FIXME: jittery)
				updated = 1;
			} // for S
			if (updated) {had_coll = 1;}
			else if (point_in_polygon_2d(pos2.x, pos2.y, &points.front(), num_sides, 0, 1)) { // test top plane (sphere on top of polygon?)
				pos2.z = p_last2.z; // assume falling/z coll
				had_coll = 1;
			}
		}
		else if (sphere_cube_intersect(pos2, radius, (*i + xlate), p_last2, p_int, cnorm, cdir, 1, xy_only)) {
			pos2 = p_int; // update current pos
			had_coll = 1; // flag as colliding, continue to look for more collisions (inside corners)
		}
	} // for i
	if (!had_coll) return 0;
	if (is_rotated()) {do_xy_rotate(rot_sin, rot_cos, center, pos2);} // rotate back
	pos = pos2;
	return had_coll;
}

unsigned building_t::check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t) const {

	if (!check_line_clip(p1-xlate, p2-xlate, bcube.d)) return 0; // no intersection
	point p1r(p1), p2r(p2);
	float tmin(0.0), tmax(1.0);
	unsigned coll(0); // 0=none, 1=side, 2=roof

	if (is_rotated()) {
		point const center(bcube.get_cube_center() + xlate);
		do_xy_rotate(-rot_sin, rot_cos, center, p1r); // inverse rotate - negate the sine term
		do_xy_rotate(-rot_sin, rot_cos, center, p2r);
	}
	p1r -= xlate; p2r -= xlate;
	float const pzmin(min(p1r.z, p2r.z)), pzmax(max(p1r.z, p2r.z));
	bool const vert(p1r.x == p2r.x && p1r.y == p2r.y);

	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (pzmin > i->d[2][1] || pzmax < i->d[2][0]) continue; // no overlap in z
		bool hit(0);

		if (use_cylinder_coll()) {
			point const cc(i->get_cube_center());
			float const dist(pt_line_dist(cc, p1r, p2r));
			vector3d const csz(i->get_size());
			float const radius(0.5*max(csz.x, csz.y));
			if (dist > radius) continue; // test conservative bounding circle
			
			if (vert) { // vertical cylinder optimization + handling of ellipsoids
				float const dx(cc.x - p1r.x), dy(cc.y - p1r.y), rx(0.5*csz.x), ry(0.5*csz.y);
				if (dx*dx/(rx*rx) + dy*dy/(ry*ry) > 1.0) continue; // no intersection (below test should return true as well)
				tmin = (i->d[2][1] - p1r.z)/(p2r.z - p1r.z);
				if (tmin < t) {t = tmin; hit = 1;}
			}
			else {
				point const cp1(cc - vector3d(0.0, 0.0, 0.5*csz.z)), cp2(cc + vector3d(0.0, 0.0, 0.5*csz.z));
				if (line_int_cylinder(p1r, p2r, cp1, cp2, radius, radius, 1, tmin) && tmin < t) {t = tmin; hit = 1;}
			}
		}
		else if (num_sides != 4) {
			vector<point> points;
			building_draw.calc_poly_pts(*this, (*i + xlate), points, num_sides);
			float const tz((i->d[2][1] - p1r.z)/(p2r.z - p1r.z)); // t value at zval = top of cube
			float const xval(p1r.x + tz*(p2r.x - p1r.x)), yval(p1r.y + tz*(p2r.y - p1r.y));

			if (point_in_polygon_2d(xval, yval, &points.front(), num_sides, 0, 1)) { // XY plane test for vertical lines and top surface
				tmin = (i->d[2][1] - p1r.z)/(p2r.z - p1r.z);
				if (tmin < t) {t = tmin; hit = 1;}
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
		else if (get_line_clip(p1r, p2r, i->d, tmin, tmax) && tmin < t) {t = tmin; hit = 1;}

		if (hit) {
			float const zval(p1.z + t*(p2.z - p1.z));
			coll = ((fabs(zval - i->d[2][1]) < 0.0001*i->get_dz()) ? 2 : 1); // test if clipped zval is close to the roof zval
		}
	} // for i
	return coll;
}

void building_t::gen_geometry(unsigned ix) {

	if (!is_valid()) return; // invalid building
	cube_t const base(parts.empty() ? bcube : parts.back());
	parts.clear(); // just in case
	building_mat_t const &mat(get_material());
	rand_gen_t rgen;
	rgen.set_state(123+ix, 345*ix);

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
	}
	if ((num_sides == 3 || num_sides == 4 || num_sides == 6) && mat.max_asf > 0.0 && rgen.rand_probability(mat.asf_prob)) { // triangles/cubes/hexagons
		alt_step_factor = max(0.0f, min(0.99f, rgen.rand_uniform(mat.min_asf, mat.max_asf)));
		if (alt_step_factor > 0.0 && !(num_sides&1)) {half_offset = 1;} // chamfered cube/hexagon
		if (alt_step_factor > 0.0) {num_sides *= 2;}
	}

	// determine the number of levels and splits
	unsigned num_levels(mat.min_levels);
	if (mat.min_levels < mat.max_levels && was_cube) {num_levels += rgen.rand()%(mat.max_levels - mat.min_levels + 1);} // only cubes are multilevel (unless min_level > 1)
	if (mat.min_level_height > 0.0) {num_levels = max(mat.min_levels, min(num_levels, unsigned(bcube.get_size().z/mat.min_level_height)));}
	num_levels = max(num_levels, 1U); // min_levels can be zero to apply more weight to 1 level buildings
	bool const do_split(num_levels < 4 && is_cube() && rgen.rand_probability(mat.split_prob)); // don't split buildings with 4 or more levels, or non-cubes

	if (num_levels == 1) { // single level
		if (do_split) {split_in_xy(base, rgen);} // generate L, T, or U shape
		else {parts.push_back(base);} // single part, entire cube
		return; // for now the bounding cube
	}
	// generate building levels and splits
	parts.resize(num_levels);
	float const height(base.get_dz()), dz(height/num_levels);

	if ((rgen.rand()&1) && !do_split) {
		point const llc(base.get_llc()), sz(base.get_size());

		for (unsigned i = 0; i < num_levels; ++i) { // generate overlapping cube levels
			cube_t &bc(parts[i]);
			bc.d[2][0] = base.d[2][0]; // z1
			bc.d[2][1] = base.d[2][0] + (i+1)*dz; // z2
			if (i > 0) {bc.d[2][1] += dz*rgen.rand_uniform(-0.5, 0.5); bc.d[2][1] = min(bc.d[2][1], base.d[2][1]);}

			for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to generate a cube that doesn't contain any existing cubes (can occasionally still fail)
				for (unsigned d = 0; d < 2; ++d) { // x,y
					bc.d[d][0] = base.d[d][0] + max(rgen.rand_uniform(-0.2, 0.45), 0.0f)*sz[d];
					bc.d[d][1] = base.d[d][1] - max(rgen.rand_uniform(-0.2, 0.45), 0.0f)*sz[d];
				}
				assert(bc.is_strictly_normalized());
				bool contains(0);
				for (unsigned j = 0; j < i; ++j) {contains |= bc.contains_cube(parts[j]);}
				if (!contains) break; // success
			} // for n
		} // for i
		return;
	}
	for (unsigned i = 0; i < num_levels; ++i) {
		cube_t &bc(parts[i]);
		if (i == 0) {bc = base;} // use full building footprint
		else {
			cube_t const &prev(parts[i-1]);
			for (unsigned d = 0; d < 2; ++d) {
				float const len(prev.d[d][1] - prev.d[d][0]);
				for (unsigned e = 0; e < 2; ++e) {
					float delta(0.0);
					if (rgen.rand()&3) {delta = rgen.rand_uniform(0.1, 0.4);} // 25% chance of no shift, 75% chance of 20-40% shift
					bc.d[d][e] = prev.d[d][e] + (e ? -delta : delta)*len;
				}
			}
			bc.d[2][0] = prev.d[2][1]; // z1
		}
		bc.d[2][1] = bc.d[2][0] + dz; // z2
		bc.normalize(); // handle XY inversion due to shift
	} // for i
	for (unsigned i = 1; i < num_levels; ++i) {
		float const ddz(rgen.rand_uniform(-0.35*dz, 0.35*dz)); // random shift in z height
		parts[i  ].d[2][0] += ddz;
		parts[i-1].d[2][1] += ddz;
	}
	if (do_split) { // generate L, T, or U shape
		cube_t const split_cube(parts.back());
		parts.pop_back();
		split_in_xy(split_cube, rgen);
	}
}

bool check_tile_smap(bool shadow_only) {
	return (!shadow_only && world_mode == WMODE_INF_TERRAIN && shadow_map_enabled());
}

void building_t::draw(shader_t &s, bool shadow_only, float far_clip, vector3d const &xlate, building_draw_t &bdraw, unsigned draw_ix) const {

	if (draw_ix == cur_draw_ix) return; // already drawn this pass
	if (!is_valid()) return; // invalid building
	cur_draw_ix = draw_ix;
	point const center(bcube.get_cube_center()), pos(center + xlate), camera(get_camera_pos());
	float const dmax(far_clip + 0.5*bcube.get_size().get_max_val());
	if (!shadow_only && !dist_less_than(camera, pos, dmax)) return; // dist clipping
	if (!camera_pdu.sphere_visible_test(pos, bcube.get_bsphere_radius())) return; // VFC
	building_mat_t const &mat(get_material());
	bool const is_close(dist_less_than(camera, pos, 0.1*far_clip));
	bool const immediate_mode(check_tile_smap(shadow_only) && try_bind_tile_smap_at_point(pos, s)); // for nearby TT tile shadow maps
	if (immediate_mode) {bdraw.begin_immediate_building();}
	vector3d view_dir(zero_vector);

	for (auto i = parts.begin(); i != parts.end(); ++i) { // multiple cubes/parts/levels case
		if (!shadow_only) {
			point ccenter(i->get_cube_center());
			if (is_rotated()) {do_xy_rotate(rot_sin, rot_cos, center, ccenter);}
			view_dir = (ccenter + xlate - camera);
		}
		bdraw.add_section(*this, *i, xlate, bcube, mat.side_tex, side_color, shadow_only, view_dir, 3); // XY
		
		if (num_sides == 4 && i->d[2][0] > bcube.d[2][0] && camera.z < i->d[2][1]) { // stacked cubes viewed from below; cur corners can have overhangs
			continue; // top surface not visible, bottom surface occluded, skip (even for shadow pass)
		}
		bdraw.add_section(*this, *i, xlate, bcube, mat.roof_tex, roof_color, shadow_only, view_dir, 4); // only Z dim
		if (is_close) {} // placeholder for drawing of building interiors, windows, detail, etc.
	}
	if (immediate_mode) {bdraw.end_immediate_building(shadow_only);}
}


unsigned const grid_sz = 32;

class building_creator_t {

	vector3d range_sz, range_sz_inv, max_extent;
	cube_t range;
	rand_gen_t rgen;
	vector<building_t> buildings;

	struct grid_elem_t {
		vector<unsigned> ixs;
		cube_t bcube;
		void add(cube_t const &c, unsigned ix) {
			if (ixs.empty()) {bcube = c;} else {bcube.union_with_cube(c);}
			ixs.push_back(ix);
		}
	};
	vector<grid_elem_t> grid;

	grid_elem_t &get_grid_elem(unsigned gx, unsigned gy) {
		assert(gx < grid_sz && gy < grid_sz);
		return grid[gy*grid_sz + gx];
	}
	grid_elem_t const &get_grid_elem(unsigned gx, unsigned gy) const {
		assert(gx < grid_sz && gy < grid_sz);
		return grid[gy*grid_sz + gx];
	}
	void get_grid_pos(point pos, unsigned ixp[2]) const { // {x,y}
		range.clamp_pt(pos);
		for (unsigned d = 0; d < 2; ++d) {
			float const v((pos[d] - range.d[d][0])*range_sz_inv[d]);
			ixp[d] = unsigned(v*(grid_sz-1));
			assert(ixp[d] < grid_sz);
		}
	}
	void get_grid_range(cube_t const &bcube, unsigned ixr[2][2]) const { // {lo,hi}x{x,y}
		get_grid_pos(bcube.get_llc(), ixr[0]);
		get_grid_pos(bcube.get_urc(), ixr[1]);
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
	vector3d const get_query_xlate() const {
		return vector3d((world_mode == WMODE_INF_TERRAIN) ? vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0) : zero_vector);
	}

public:
	building_creator_t() : max_extent(zero_vector) {}
	bool empty() const {return buildings.empty();}
	void clear() {buildings.clear(); grid.clear();}
	vector3d const &get_max_extent() const {return max_extent;}
	building_t const &get_building(unsigned ix) const {assert(ix < buildings.size()); return buildings[ix];}

	void gen(building_params_t const &params) {
		if (params.materials.empty()) return; // no materials
		timer_t timer("Gen Buildings");
		float const def_water_level(get_water_z_height());
		vector3d const offset(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0);
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? offset : zero_vector); // cancel out xoff2/yoff2 translate
		range = params.materials.front().pos_range; // range is union over all material ranges
		for (auto i = params.materials.begin()+1; i != params.materials.end(); ++i) {range.union_with_cube(i->pos_range);}
		range     += ((world_mode == WMODE_INF_TERRAIN) ? zero_vector : offset);
		range_sz   = range.get_size(); // Note: place_radius is relative to range cube center
		max_extent = zero_vector;
		UNROLL_3X(range_sz_inv[i_] = 1.0/range_sz[i_];)
		clear();
		buildings.reserve(params.num_place);
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_tries(0), num_gen(0), num_skip(0);
		rgen.set_state(rand_gen_index, 123); // update when mesh changes, otherwise determinstic

		for (unsigned i = 0; i < params.num_place; ++i) {
			point center(all_zeros);
			
			for (unsigned n = 0; n < params.num_tries; ++n) { // 10 tries to find a non-overlapping building placement
				building_t b(params.choose_rand_mat(rgen)); // set material
				building_mat_t const &mat(b.get_material());
				point const place_center(mat.pos_range.get_cube_center());
				bool keep(0);

				for (unsigned m = 0; m < params.num_tries; ++m) {
					for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(mat.pos_range.d[d][0], mat.pos_range.d[d][1]);} // x,y
					if (mat.place_radius == 0.0 || dist_xy_less_than(center, place_center, mat.place_radius)) {keep = 1; break;}
				}
				if (!keep) continue; // placement failed, skip
				center.z = get_exact_zval(center.x+xlate.x, center.y+xlate.y);

				for (unsigned d = 0; d < 3; ++d) { // x,y,z
					float const sz(0.5*rgen.rand_uniform(mat.sz_range.d[d][0], mat.sz_range.d[d][1]));
					b.bcube.d[d][0] = center[d] - ((d == 2) ? 0.0 : sz); // only in XY
					b.bcube.d[d][1] = center[d] + sz;
				} // for d
				++num_tries;
				float const z_sea_level(center.z - def_water_level);
				if (z_sea_level < 0.0) break; // skip underwater buildings, failed placement
				if (z_sea_level < mat.min_alt || z_sea_level > mat.max_alt) break;
				b.gen_rotation(rgen);
				++num_gen;

				// check building for overlap with other buildings
				float const expand(b.is_rotated() ? 0.05 : 0.1); // expand by 5-10%
				cube_t test_bc(b.bcube);
				test_bc.expand_by(expand*b.bcube.get_size());
				sphere_t const bsphere(test_bc.get_bcylin()); // only care about XY radius
				bool overlaps(0);
				unsigned ixr[2][2];
				get_grid_range(b.bcube, ixr);

				for (unsigned y = ixr[0][1]; y <= ixr[1][1] && !overlaps; ++y) {
					for (unsigned x = ixr[0][0]; x <= ixr[1][0] && !overlaps; ++x) {
						grid_elem_t const &ge(get_grid_elem(x, y));
						if (!test_bc.intersects_xy(ge.bcube)) continue;

						for (auto g = ge.ixs.begin(); g != ge.ixs.end(); ++g) {
							building_t const &ob(get_building(*g));
							if (test_bc.intersects_xy(ob.bcube) && ob.check_bcube_overlap_xy(b, expand)) {overlaps = 1; break;}
						}
					} // for x
				} // for y
				if (!overlaps) {
					mat.side_color.gen_color(b.side_color, rgen);
					mat.roof_color.gen_color(b.roof_color, rgen);
					add_to_grid(b.bcube, buildings.size());
					vector3d const sz(b.bcube.get_size());
					float const mult[3] = {0.5, 0.5, 1.0}; // half in X,Y and full in Z
					UNROLL_3X(max_extent[i_] = max(max_extent[i_], mult[i_]*sz[i_]);)
					buildings.push_back(b);
					break; // done
				}
			} // for n
		} // for i
		timer.end();

		if (params.flatten_mesh) {
			timer_t timer("Gen Building Zvals");
			bool const do_flatten(using_tiled_terrain_hmap_tex());

#pragma omp parallel for schedule(static,1)
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
					else if (!b.parts.empty()) {b.parts.back().d[2][0] = b.bcube.d[2][0];} // update base z1
				}
			} // for i
			if (do_flatten) { // use conservative zmin for grid
				for (auto i = grid.begin(); i != grid.end(); ++i) {i->bcube.d[2][0] = def_water_level;}
			}
		} // if flatten_mesh

		timer_t timer2("Gen Building Geometry");
#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < (int)buildings.size(); ++i) {buildings[i].gen_geometry(i);}

		cout << "WM: " << world_mode << " Buildings: " << params.num_place << " / " << num_tries << " / " << num_gen
			 << " / " << buildings.size() << " / " << (buildings.size() - num_skip) << endl;
	}

	void draw(bool shadow_only, vector3d const &xlate) const {
		if (empty()) return;
		//timer_t timer("Draw Buildings"); // 1.7ms, 2.3ms with shadow maps, 2.8ms with AO, 3.3s with rotations (currently 2.5)
		float const far_clip(get_inf_terrain_fog_dist());
		point const camera(get_camera_pos());
		int const use_bmap(global_building_params.has_normal_map);
		bool const use_tt_smap(check_tile_smap(shadow_only));
		static unsigned draw_ix(0); ++draw_ix;
		shader_t s;
		fgPushMatrix();
		translate_to(xlate);
		if (!shadow_only) {building_draw.init_draw_frame();}

		if (use_tt_smap) { // pre-pass to render buildings in nearby tiles that have shadow maps
			setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 1, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
			s.add_uniform_float("z_bias", cobj_z_bias);
			s.add_uniform_float("pcf_offset", 10.0*shadow_map_pcf_offset);
		}
		for (auto g = grid.begin(); g != grid.end(); ++g) {
			point const pos(g->bcube.get_cube_center() + xlate);
			if (!shadow_only && !dist_less_than(camera, pos, (far_clip + 0.5*g->bcube.get_size().get_max_val()))) continue; // too far
			if (!camera_pdu.sphere_visible_test(pos, g->bcube.get_bsphere_radius())) continue; // VFC
			for (auto i = g->ixs.begin(); i != g->ixs.end(); ++i) {buildings[*i].draw(s, shadow_only, far_clip, xlate, building_draw, draw_ix);}
		}
		if (use_tt_smap) {s.end_shader();}
		if (shadow_only) {s.begin_color_only_shader();} // really don't even need colors
		else { // main/batched draw pass
			bool const v(world_mode == WMODE_GROUND), indir(v), dlights(v), use_smap(v);
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, dlights, 0, 0, use_smap, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
		}
		building_draw.draw_and_clear(shadow_only);
		s.end_shader();
		fgPopMatrix();
	}

	bool check_sphere_coll(point &pos, point const &p_last, float radius, bool xy_only=0) const {
		if (empty()) return 0;
		vector3d const xlate(get_query_xlate());
		cube_t bcube; bcube.set_from_sphere((pos - xlate), radius);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		float const dist(p2p_dist(pos, p_last));

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (!sphere_cube_intersect(pos, (radius + dist), (ge.bcube + xlate))) continue;

				// Note: assumes buildings are separated so that only one sphere collision can occur
				for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) {
					if (get_building(*b).check_sphere_coll(pos, p_last, xlate, radius, xy_only)) return 1;
				} // for g
			} // for x
		} // for y
		return 0;
	}

	unsigned check_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix) const {
		if (empty()) return 0;
		bool const vertical(p1.x == p2.x && p1.y == p2.y);
		vector3d const xlate(get_query_xlate());
		cube_t bcube(p1-xlate, p2-xlate);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		point end_pos(p2);
		unsigned coll(0); // 0=none, 1=side, 2=roof
		t = 1.0; // start at end point

		// for now, just do a slow iteration over every grid element within the line's bbox in XY
		// Note: should probably iterate over the grid in XY order from the start to the end of the line, or better yet use a line drawing algorithm
		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (!check_line_clip(p1, end_pos, ge.bcube.d)) continue; // no intersection - skip this grid

				for (auto b = ge.ixs.begin(); b != ge.ixs.end(); ++b) { // Note: okay to check the same building more than once
					building_t const &building(get_building(*b));
					if (!building.bcube.intersects(bcube)) continue;
					float t_new(t);
					unsigned const ret(building.check_line_coll(p1, p2, xlate, t_new));

					if (ret && t_new <= t) { // closer hit pos, update state
						t       = t_new;
						hit_bix = *b;
						coll    = ret;
						end_pos = p1 + t*(p2 - p1);
						if (vertical) return coll; // vertical lines can only intersect one building
					}
				}
			} // for x
		} // for y
		return coll;
	}
};


building_creator_t building_creator;

void gen_buildings() {building_creator.gen(global_building_params);}
void draw_buildings(bool shadow_only, vector3d const &xlate) {building_creator.draw(shadow_only, xlate);}

bool check_buildings_point_coll(point const &pos, bool apply_tt_xlate, bool xy_only) {
	return check_buildings_sphere_coll(pos, 0.0, apply_tt_xlate, xy_only);
}
bool check_buildings_sphere_coll(point const &pos, float radius, bool apply_tt_xlate, bool xy_only) {
	point center(pos);
	if (apply_tt_xlate && world_mode == WMODE_INF_TERRAIN) {center += vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0);} // apply xlate for all static objects - not the camera
	return building_creator.check_sphere_coll(center, pos, radius, xy_only);
}
bool proc_buildings_sphere_coll(point &pos, point const &p_int, float radius, bool xy_only) {
	return building_creator.check_sphere_coll(pos, p_int, radius, xy_only);
}
unsigned check_buildings_line_coll(point const &p1, point const &p2, float &t, unsigned &hit_bix, bool apply_tt_xlate) {
	vector3d const xlate((apply_tt_xlate && world_mode == WMODE_INF_TERRAIN) ? vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0) : zero_vector);
	return building_creator.check_line_coll(p1+xlate, p2+xlate, t, hit_bix);
}
bool get_buildings_line_hit_color(point const &p1, point const &p2, colorRGBA &color) {
	float t(0.0); // unused
	unsigned hit_bix(0);
	unsigned const ret(check_buildings_line_coll(p1, p2, t, hit_bix, 1)); // apply_tt_xlate=1; 0=no hit, 1=hit side, 2=hit roof
	if (ret == 0) return 0;
	building_t const &b(building_creator.get_building(hit_bix));
	color = ((ret == 2) ? b.get_avg_roof_color() : b.get_avg_side_color());
	return 1;
}
vector3d const &get_buildings_max_extent() {return building_creator.get_max_extent();} // used for TT shadow bounds


