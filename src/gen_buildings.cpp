// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "file_utils.h"

using std::string;

extern float water_plane_z;

// TODO:
// multiple texture/color sets for sides/roof


struct tid_nm_pair_t {

	int tid, nm_tid;
	float tscale;

	tid_nm_pair_t() : tid(-1), nm_tid(-1), tscale(1.0) {}
	bool enabled() const {return (tid >= 0 || nm_tid >= 0);}
	bool operator==(tid_nm_pair_t const &t) const {return (tid == t.tid && nm_tid == t.nm_tid && tscale == t.tscale);}

	void set_gl() const {
		select_texture(tid);
		select_multitex(nm_tid, 5);
	}
};

struct building_mat_t {

	tid_nm_pair_t side_tex, roof_tex;
};

struct building_params_t : public building_mat_t {

	bool flatten_mesh;
	unsigned num;
	cube_t sz_range, pos_range; // z is unused?
	colorRGBA color_min, color_max; // alpha is unused?

	building_params_t(unsigned num_=0) : flatten_mesh(0), num(num_), sz_range(all_zeros), pos_range(all_zeros), color_min(WHITE), color_max(WHITE) {}
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

	if (str == "flatten_mesh") {
		if (!read_bool(fp, global_building_params.flatten_mesh)) {buildings_file_err("flatten_mesh", error);}
	}
	else if (str == "num") {
		if (!read_uint(fp, global_building_params.num)) {buildings_file_err("num", error);}
	}
	else if (str == "side_tscale") {
		if (!read_float(fp, global_building_params.side_tex.tscale)) {buildings_file_err("side_tscale", error);}
	}
	else if (str == "roof_tscale") {
		if (!read_float(fp, global_building_params.roof_tex.tscale)) {buildings_file_err("roof_tscale", error);}
	}
	else if (str == "side_tid") {
		if (!read_str(fp, strc)) {buildings_file_err("side_tid", error);}
		global_building_params.side_tex.tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "side_nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err("side_nm_tid", error);}
		global_building_params.side_tex.nm_tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "roof_tid") {
		if (!read_str(fp, strc)) {buildings_file_err("roof_tid", error);}
		global_building_params.roof_tex.tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "roof_nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err("roof_nm_tid", error);}
		global_building_params.roof_tex.nm_tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "size_range") {
		if (!read_cube(fp, global_building_params.sz_range)) {buildings_file_err("size_range", error);}
	}
	else if (str == "pos_range") {
		if (!read_cube(fp, global_building_params.pos_range)) {buildings_file_err("pos_range", error);}
	}
	else if (str == "color_min") {
		if (!read_color(fp, global_building_params.color_min)) {buildings_file_err("color_min", error);}
	}
	else if (str == "color_max") {
		if (!read_color(fp, global_building_params.color_max)) {buildings_file_err("color_max", error);}
	}
	else {
		cout << "Unrecognized buildings keyword in input file: " << str << endl;
		error = 1;
	}
	return !error;
}


struct building_t : public building_mat_t {

	colorRGBA color;
	cube_t bcube;

	building_t(building_mat_t const &mat=building_mat_t()) : building_mat_t(mat), color(WHITE) {}

	void draw(shader_t &s, bool shadow_only, float far_clip, vector3d const &xlate=zero_vector) const {
		// FIXME: more efficient, use batched verts or VBO, cache verts, store color in verts
		point const center(bcube.get_cube_center()), pos(center + xlate);
		float const dmax(far_clip + 0.5*bcube.get_size().get_max_val());
		if (!dist_less_than(get_camera_pos(), pos, dmax)) return; // dist clipping
		if (!camera_pdu.sphere_visible_test(pos, bcube.get_bsphere_radius())) return; // VFC
		vector3d const sz(bcube.get_size()), view_dir(pos - get_camera_pos());
		bool const single_pass(shadow_only || side_tex == roof_tex);
		if (!shadow_only) {s.set_cur_color(color);}
		
		// draw sides
		if (!shadow_only) {side_tex.set_gl();}
		draw_cube(center, sz.x, sz.y, sz.z, (!shadow_only && side_tex.enabled()), 0, side_tex.tscale, 1, &view_dir, (single_pass ? 7 : 3), 1);
		
		if (!single_pass) { // draw roof (and floor if at water edge)
			roof_tex.set_gl();
			draw_cube(center, sz.x, sz.y, sz.z, roof_tex.enabled(), 0, roof_tex.tscale, 1, &view_dir, 4, 1); // only Z dim
		}
	}
};

unsigned const grid_sz = 32;

class building_creator_t {

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
	bool empty() const {return buildings.empty();}
	void clear() {buildings.clear(); grid.clear();}

	void gen(building_params_t const &params) {
		timer_t timer("Gen Buildings");
		vector3d const xlate((world_mode == WMODE_INF_TERRAIN) ? vector3d(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0) : zero_vector); // cancel out xoff2/yoff2 translate
		range    = params.pos_range;
		range_sz = range.get_size();
		UNROLL_3X(range_sz_inv[i_] = 1.0/range_sz[i_];)
		clear();
		buildings.reserve(params.num);
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_tries(0), num_gen(0);

		for (unsigned i = 0; i < params.num; ++i) {
			building_t b(params); // copy material
			point center;
			
			for (unsigned n = 0; n < 10; ++n) { // 10 tries to find a non-overlapping building placement
				for (unsigned d = 0; d < 3; ++d) { // x,y,z
					if (d < 2) {center[d] = rgen.rand_uniform(range.d[d][0], range.d[d][1]);} // x,y
					else {center[d] = get_exact_zval(center.x+xlate.x, center.y+xlate.y);} // z
					float const sz(0.5*rgen.rand_uniform(params.sz_range.d[d][0], params.sz_range.d[d][1]));
					b.bcube.d[d][0] = center[d] - ((d == 2) ? 0.0 : sz); // only in XY
					b.bcube.d[d][1] = center[d] + sz;
				} // for d
				++num_tries;
				if (center.z < water_plane_z) break; // skip underwater buildings, failed placement
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
					for (unsigned d = 0; d < 4; ++d) {b.color[d] = rgen.rand_uniform(params.color_min[d], params.color_max[d]);}
					add_to_grid(b.bcube, buildings.size());
					buildings.push_back(b);
					break; // done
				}
			} // for n
		} // for i
		timer.end();
		cout << "Buildings: " << params.num << " / " << num_tries << " / " << num_gen << " / " << buildings.size() << endl;

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
					for (int d = 0; d < 4; ++d) {zmin = min(zmin, get_exact_zval(b.bcube.d[0][d&1]+xlate.x, b.bcube.d[1][d>>1]+xlate.y));}
					zmin = max(zmin, water_plane_z); // don't go below the water
				}
			} // for i
			if (do_flatten) { // use conservative zmin for grid
				for (auto i = grid.begin(); i != grid.end(); ++i) {i->bcube.d[2][0] = water_plane_z;}
			}
		}
	}

	void draw(bool shadow_only, vector3d const &xlate) const {
		if (empty()) return;
		//timer_t timer("Draw Buildings");
		fgPushMatrix();
		translate_to(xlate);
		shader_t s;

		if (shadow_only) {
			s.begin_color_only_shader(); // really don't even need colors
		}
		else {
			int const use_bmap(global_building_params.side_tex.nm_tid >= 0 || global_building_params.roof_tex.nm_tid >= 0), is_outside(1);
			bool const indir(0), use_smap(0); // FIXME: shadows between buildings
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, 0, 0, 0, use_smap, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, is_outside);
		}
		float const far_clip(get_inf_terrain_fog_dist());
		for (auto i = buildings.begin(); i != buildings.end(); ++i) {i->draw(s, shadow_only, far_clip, xlate);}
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
					point p_int;
					vector3d cnorm; // unused
					unsigned cdir(0); // unused

					if (sphere_cube_intersect(pos, radius, (buildings[*g].bcube + xlate), p_last, p_int, cnorm, cdir, 1, 0)) {
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


