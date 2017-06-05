// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "file_utils.h"

using std::string;

extern int rand_gen_index, display_mode;

// TODO:
// windows in brick/block buildings
// TT tile shadows
// L-shaped/non-rectangular buildings

struct tid_nm_pair_t {

	int tid, nm_tid; // Note: assumes each tid has only one nm_tid
	float tscale; // FIXME: tscale_x vs. tscale_y?

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale(1.0) {}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}
	bool operator==(tid_nm_pair_t const &t) const {return (tid == t.tid && nm_tid == t.nm_tid && tscale == t.tscale);}

	void set_gl() const {
		select_texture(tid);
		select_multitex(nm_tid, 5);
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

	unsigned min_levels, max_levels;
	color_range_t side_color, roof_color;
	cube_t sz_range;

	building_mat_t() : min_levels(1), max_levels(1), sz_range(1,1,1,1,1,1) {}
	bool has_normal_map() const {return (side_tex.nm_tid >= 0 || roof_tex.nm_tid >= 0);}
};

struct building_params_t {

	bool flatten_mesh, has_normal_map;
	unsigned num_place, num_tries, cur_prob;
	float place_radius, max_delta_z;
	cube_t pos_range; // z is unused?
	building_mat_t cur_mat;
	vector<building_mat_t> materials;
	vector<unsigned> mat_gen_ix;

	building_params_t(unsigned num_place_=0) : flatten_mesh(0), has_normal_map(0), num_place(num_place_), num_tries(10),
		cur_prob(1), place_radius(0.0), max_delta_z(0.0), pos_range(-100,100,-100,100,0,0) {}
	
	void add_cur_mat() {
		unsigned const mat_ix(materials.size());
		for (unsigned n = 0; n < cur_prob; ++n) {mat_gen_ix.push_back(mat_ix);} // add more references to this mat for higher probability
		materials.push_back(cur_mat);
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
	else if (str == "pos_range") {
		if (!read_cube(fp, global_building_params.pos_range)) {buildings_file_err(str, error);}
	}
	else if (str == "place_radius") {
		if (!read_float(fp, global_building_params.place_radius)) {buildings_file_err(str, error);}
	}
	else if (str == "max_delta_z") {
		if (!read_float(fp, global_building_params.max_delta_z)) {buildings_file_err(str, error);}
	}
	// material parameters
	else if (str == "min_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.min_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "max_levels") {
		if (!read_uint(fp, global_building_params.cur_mat.max_levels)) {buildings_file_err(str, error);}
	}
	else if (str == "size_range") {
		if (!read_cube(fp, global_building_params.cur_mat.sz_range)) {buildings_file_err(str, error);}
	}
	// material textures
	else if (str == "side_tscale") {
		if (!read_float(fp, global_building_params.cur_mat.side_tex.tscale)) {buildings_file_err(str, error);}
	}
	else if (str == "roof_tscale") {
		if (!read_float(fp, global_building_params.cur_mat.roof_tex.tscale)) {buildings_file_err(str, error);}
	}
	else if (str == "side_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_tex.tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "side_nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.side_tex.nm_tid = get_texture_by_name(std::string(strc), 1);
	}
	else if (str == "roof_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_tex.tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "roof_nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err(str, error);}
		global_building_params.cur_mat.roof_tex.nm_tid = get_texture_by_name(std::string(strc), 1);
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

struct building_t {

	unsigned mat_ix;
	colorRGBA side_color, roof_color;
	cube_t bcube;
	vector<cube_t> levels;
	mutable unsigned cur_draw_ix;

	building_t(unsigned mat_ix_=0) : mat_ix(mat_ix_), side_color(WHITE), roof_color(WHITE), cur_draw_ix(0) {bcube.set_to_zeros();}
	bool is_valid() const {return !bcube.is_all_zeros();}
	building_mat_t const &get_material() const {return global_building_params.get_material(mat_ix);}
	void gen_levels(unsigned ix);
	void draw(bool shadow_only, float far_clip, vector3d const &xlate, building_draw_t &bdraw, unsigned draw_ix) const;
};


class building_draw_t {

	struct draw_block_t {
		tid_nm_pair_t tex;
		vector<vert_norm_comp_tc_color> verts;

		void draw_and_clear(bool shadow_only) {
			if (verts.empty()) return;
			if (!shadow_only) {tex.set_gl();}
			draw_quad_verts_as_tris(verts);
			verts.clear();
		}
		bool empty() const {return verts.empty();}
	};
	vector<draw_block_t> to_draw; // one per texture, assumes tids are dense

	vector<vert_norm_comp_tc_color> &get_verts(tid_nm_pair_t const &tex) {
		unsigned const ix((tex.tid >= 0) ? (tex.tid+1) : 0);
		if (ix >= to_draw.size()) {to_draw.resize(ix+1);}
		if (to_draw[ix].empty()) {to_draw[ix].tex = tex;} // copy material first time
		else {assert(to_draw[ix].tex.nm_tid == tex.nm_tid);} // else normal maps must agree
		return to_draw[ix].verts;
	}
	void add_cube_verts(cube_t const &cube, vector<vert_norm_comp_tc_color> &verts, colorRGBA const &color,
		bool texture, float texture_scale, vector3d const *const view_dir, unsigned dim_mask)
	{
		texture_scale *= 2.0; // adjust for local vs. global space change
		point const scale(cube.get_size());
		vector3d const xlate(cube.get_llc()); // move origin from center to min corner
		vert_norm_comp_tc_color vert;
		vert.set_c4(color); // color is shared across all verts

		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);
			if (!(dim_mask & (1<<n))) continue;

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (view_dir && (((*view_dir)[n] < 0.0) ^ j)) continue; // back facing
				vert.n[d[0]] = 0;
				vert.n[d[1]] = 0;
				vert.n[n] = (j ? 127 : -128); // -1.0 or 1.0
				point pt;
				pt[n] = j;

				for (unsigned s1 = 0; s1 < 2; ++s1) {
					pt[d[1]] = s1;

					for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
						pt[d[0]] = k^j^s1^1; // need to orient the vertices differently for each side
						vert.v   = pt*scale + xlate;

						if (texture) {
							bool const st(i&1);
							vert.t[ st] = texture_scale*vert.v[d[1]];
							vert.t[!st] = texture_scale*vert.v[d[0]];
						}
						verts.push_back(vert);
					}
				}
			} // for j
		} // for i
	}
public:
	void add_cube(cube_t const &cube, tid_nm_pair_t const &tex, colorRGBA const &color, bool shadow_only, vector3d const *const view_dir, unsigned dim_mask) {
		add_cube_verts(cube, get_verts(tex), color, (!shadow_only && tex.enabled()), tex.tscale, view_dir, dim_mask);
	}
	void draw_and_clear(bool shadow_only) {
		for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->draw_and_clear(shadow_only);}
	}
};

building_draw_t building_draw;


void building_t::gen_levels(unsigned ix) {

	if (!is_valid()) return; // invalid building
	building_mat_t const &mat(get_material());
	// use ix value as the seed/hash; at least one level
	unsigned num_levels(mat.min_levels);
	if (mat.min_levels < mat.max_levels) {num_levels += ix%(mat.max_levels - mat.min_levels + 1);}
	num_levels = max(num_levels, 1U); // min_levels can be zero to apply more weight to 1 level buildings
	if (num_levels == 1) return; // single level, for now the bounding cube
	levels.resize(num_levels);
	rand_gen_t rgen;
	rgen.set_state(ix, 345);
	float const height(bcube.d[2][1] - bcube.d[2][0]), dz(height/num_levels);

	for (unsigned i = 0; i < num_levels; ++i) {
		cube_t &bc(levels[i]);
		if (i == 0) {bc = bcube;} // use full building footprint
		else {
			cube_t const &prev(levels[i-1]);
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
		levels[i  ].d[2][0] += ddz;
		levels[i-1].d[2][1] += ddz;
	}
}

void building_t::draw(bool shadow_only, float far_clip, vector3d const &xlate, building_draw_t &bdraw, unsigned draw_ix) const {

	// store in VBO?
	if (!is_valid()) return; // invalid building
	if (draw_ix == cur_draw_ix) return; // already drawn this pass
	cur_draw_ix = draw_ix;
	point const center(bcube.get_cube_center()), pos(center + xlate), camera(get_camera_pos());
	float const dmax(far_clip + 0.5*bcube.get_size().get_max_val());
	if (!shadow_only && !dist_less_than(camera, pos, dmax)) return; // dist clipping
	if (!camera_pdu.sphere_visible_test(pos, bcube.get_bsphere_radius())) return; // VFC
	building_mat_t const &mat(get_material());

	if (levels.empty()) { // single cube/level case
		vector3d const view_dir(pos - camera);
		vector3d const *const vdir(shadow_only ? nullptr : &view_dir);
		bdraw.add_cube(bcube, mat.side_tex, side_color, shadow_only, vdir, 3); // XY sides
		bdraw.add_cube(bcube, mat.roof_tex, roof_color, shadow_only, vdir, 4); // Z roof (and floor if at water edge)
	}
	for (auto i = levels.begin(); i != levels.end(); ++i) { // multiple cubes/levels case
		vector3d const view_dir(shadow_only ? zero_vector : (i->get_cube_center() + xlate - camera));
		vector3d const *const vdir(shadow_only ? nullptr : &view_dir);
		bdraw.add_cube(*i, mat.side_tex, side_color, shadow_only, vdir, 3); // XY
		if (i != levels.begin() && camera.z < i->d[2][1]) continue; // top surface not visible, bottom surface occluded, skip (even for shadow pass)
		bdraw.add_cube(*i, mat.roof_tex, roof_color, shadow_only, vdir, 4); // only Z dim
	}
}


unsigned const grid_sz = 32;

class building_creator_t {

	float place_radius;
	vector3d range_sz, range_sz_inv;
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

public:
	building_creator_t() : place_radius(0.0) {}
	bool empty() const {return buildings.empty();}
	void clear() {buildings.clear(); grid.clear();}

	void gen(building_params_t const &params) {
		timer_t timer("Gen Buildings");
		float const def_water_level(get_water_z_height());
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? vector3d(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0) : zero_vector); // cancel out xoff2/yoff2 translate
		range    = params.pos_range - xlate;
		range_sz = range.get_size();
		place_radius = params.place_radius;
		UNROLL_3X(range_sz_inv[i_] = 1.0/range_sz[i_];)
		clear();
		buildings.reserve(params.num_place);
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_tries(0), num_gen(0), num_skip(0);
		point const place_center(range.get_cube_center());
		rgen.set_state(rand_gen_index, 123); // update when mesh changes, otherwise determinstic

		for (unsigned i = 0; i < params.num_place; ++i) {
			building_t b(params.choose_rand_mat(rgen)); // set material
			building_mat_t const &mat(b.get_material());
			point center(all_zeros);
			
			for (unsigned n = 0; n < params.num_tries; ++n) { // 10 tries to find a non-overlapping building placement
				bool keep(0);

				for (unsigned m = 0; m < params.num_tries; ++m) {
					for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(range.d[d][0], range.d[d][1]);} // x,y
					if (place_radius == 0.0 || dist_xy_less_than(center, place_center, place_radius)) {keep = 1; break;}
				}
				if (!keep) continue; // placement failed, skip
				center.z = get_exact_zval(center.x+xlate.x, center.y+xlate.y);

				for (unsigned d = 0; d < 3; ++d) { // x,y,z
					float const sz(0.5*rgen.rand_uniform(mat.sz_range.d[d][0], mat.sz_range.d[d][1]));
					b.bcube.d[d][0] = center[d] - ((d == 2) ? 0.0 : sz); // only in XY
					b.bcube.d[d][1] = center[d] + sz;
				} // for d
				++num_tries;
				if (center.z < def_water_level) break; // skip underwater buildings, failed placement
				++num_gen;

				// check building for overlap with other buildings
				cube_t test_bc(b.bcube);
				test_bc.expand_by(0.1*b.bcube.get_size()); // expand by 10%
				bool overlaps(0);
				unsigned ixr[2][2];
				get_grid_range(b.bcube, ixr);

				for (unsigned y = ixr[0][1]; y <= ixr[1][1] && !overlaps; ++y) {
					for (unsigned x = ixr[0][0]; x <= ixr[1][0] && !overlaps; ++x) {
						grid_elem_t const &ge(get_grid_elem(x, y));
						if (!test_bc.intersects_xy(ge.bcube)) continue;

						for (auto g = ge.ixs.begin(); g != ge.ixs.end(); ++g) {
							assert(*g < buildings.size());
							if (test_bc.intersects_xy(buildings[*g].bcube)) {overlaps = 1; break;} // Note: only check for XY intersection
						}
					}
				}
				if (!overlaps) {
					mat.side_color.gen_color(b.side_color, rgen);
					mat.roof_color.gen_color(b.roof_color, rgen);
					add_to_grid(b.bcube, buildings.size());
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
					flatten_hmap_region(b.bcube); // flatten the mesh under the bcube to a height of mesh_zval
				}
				else { // extend building bottom downward to min mesh height
					float &zmin(b.bcube.d[2][0]); // Note: grid bcube z0 value won't be correct, but will be fixed conservatively below
					float const zmin0(zmin);
					unsigned num_below(0);
					
					for (int d = 0; d < 4; ++d) {
						float const zval(get_exact_zval(b.bcube.d[0][d&1]+xlate.x, b.bcube.d[1][d>>1]+xlate.y));
						zmin = min(zmin, zval);
						num_below += (zval < def_water_level);
					}
					zmin = max(zmin, def_water_level); // don't go below the water
					if (num_below > 2 || // more than 2 corners underwater
						(params.max_delta_z > 0.0 && (zmin0 - zmin) > params.max_delta_z)) // too steep of a slope
					{
						b.bcube.set_to_zeros();
						++num_skip;
					}
				}
			} // for i
			if (do_flatten) { // use conservative zmin for grid
				for (auto i = grid.begin(); i != grid.end(); ++i) {i->bcube.d[2][0] = def_water_level;}
			}
		} // if flatten_mesh

		timer_t timer2("Gen Building Levels");
#pragma omp parallel for schedule(static,1)
		for (int i = 0; i < (int)buildings.size(); ++i) {buildings[i].gen_levels(i);}

		cout << "Buildings: " << params.num_place << " / " << num_tries << " / " << num_gen << " / " << buildings.size() << " / " << (buildings.size() - num_skip) << endl;
	}

	void draw(bool shadow_only, vector3d const &xlate) const {
		if (empty()) return;
		//timer_t timer("Draw Buildings"); // 1.8ms for <=10 level buildings
		fgPushMatrix();
		translate_to(xlate);
		shader_t s;

		if (shadow_only) {
			s.begin_color_only_shader(); // really don't even need colors
		}
		else {
			int const use_bmap(global_building_params.has_normal_map), is_outside(1);
			bool const v(world_mode == WMODE_GROUND), indir(v), dlights(v), use_smap(v); // shadows between buildings in TT mode?
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, dlights, 0, 0, use_smap, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, is_outside);
		}
		float const far_clip(get_inf_terrain_fog_dist());
		point const camera(get_camera_pos());
		static unsigned draw_ix(0); ++draw_ix;

		for (auto g = grid.begin(); g != grid.end(); ++g) {
			point const pos(g->bcube.get_cube_center() + xlate);
			if (!shadow_only && !dist_less_than(camera, pos, (far_clip + 0.5*g->bcube.get_size().get_max_val()))) continue; // too far
			if (!camera_pdu.sphere_visible_test(pos, g->bcube.get_bsphere_radius())) continue; // VFC
			for (auto i = g->ixs.begin(); i != g->ixs.end(); ++i) {buildings[*i].draw(shadow_only, far_clip, xlate, building_draw, draw_ix);}
		}
		building_draw.draw_and_clear(shadow_only);
		s.end_shader();
		fgPopMatrix();
	}

	bool check_sphere_coll(point &pos, point const &p_last, float radius=0.0) const {
		if (empty()) return 0;
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0) : zero_vector);
		cube_t bcube; bcube.set_from_sphere((pos - xlate), radius);
		unsigned ixr[2][2];
		get_grid_range(bcube, ixr);
		float const dist(p2p_dist(pos, p_last));

		for (unsigned y = ixr[0][1]; y <= ixr[1][1]; ++y) {
			for (unsigned x = ixr[0][0]; x <= ixr[1][0]; ++x) {
				grid_elem_t const &ge(get_grid_elem(x, y));
				if (!sphere_cube_intersect(pos, (radius + dist), (ge.bcube + xlate))) continue;

				for (auto g = ge.ixs.begin(); g != ge.ixs.end(); ++g) {
					assert(*g < buildings.size());
					cube_t const bc(buildings[*g].bcube);
					if (bc.is_all_zeros()) continue; // invalid building
					point p_int;
					vector3d cnorm; // unused
					unsigned cdir(0); // unused

					if (sphere_cube_intersect(pos, radius, (bc + xlate), p_last, p_int, cnorm, cdir, 1, 0)) {
						pos = p_int;
						return 1; // Note: assumes buildings are separated so that only one sphere collision can occur
					}
				}
			}
		}
		return 0;
	}
};


building_creator_t building_creator;

void gen_buildings() {building_creator.gen(global_building_params);}
void draw_buildings(bool shadow_only, vector3d const &xlate) {building_creator.draw(shadow_only, xlate);}
bool check_buildings_point_coll(point const &pos) {return check_buildings_sphere_coll(pos, 0.0);}
bool check_buildings_sphere_coll(point const &pos, float radius) {point pos2(pos); return building_creator.check_sphere_coll(pos2, pos, radius);}
bool proc_buildings_sphere_coll(point &pos, point const &p_int, float radius) {return building_creator.check_sphere_coll(pos, p_int, radius);}


