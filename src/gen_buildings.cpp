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


struct building_params_t {

	bool flatten_mesh;
	unsigned num;
	int tid, nm_tid;
	cube_t sz_range, pos_range; // z is unused?
	colorRGBA color_min, color_max; // alpha is unused?

	building_params_t(unsigned num_=0) : flatten_mesh(0), num(num_), tid(-1), nm_tid(-1), sz_range(all_zeros), pos_range(all_zeros), color_min(WHITE), color_max(WHITE) {}
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
	else if (str == "tid") {
		if (!read_str(fp, strc)) {buildings_file_err("tid", error);}
		global_building_params.tid = get_texture_by_name(std::string(strc));
	}
	else if (str == "nm_tid") {
		if (!read_str(fp, strc)) {buildings_file_err("nm_tid", error);}
		global_building_params.nm_tid = get_texture_by_name(std::string(strc));
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


struct building_t {

	int tid, nm_tid;
	colorRGBA color;
	cube_t bcube;

	building_t() : tid(-1), nm_tid(-1), color(WHITE) {}

	void draw(shader_t &s, bool shadow_only, vector3d const &xlate=zero_vector) const { // FIXME: more efficient, use batched verts or VBO
		point const pos(bcube.get_cube_center() + xlate);
		if (world_mode == WMODE_INF_TERRAIN && distance_to_camera(pos) > get_inf_terrain_fog_dist() + bcube.get_size().get_max_val()) return; // dist clipping
		if (!camera_pdu.sphere_and_cube_visible_test(pos, bcube.get_bsphere_radius(), (bcube + xlate))) return; // FIXME: cache center and radius?
		if (!shadow_only) {select_texture(tid); select_multitex(nm_tid, 5);}
		s.set_cur_color(color);
		draw_simple_cube(bcube, (!shadow_only && tid >= 0));
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
		vector3d const xlate(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0); // cancel out xoff2/yoff2 translate
		range    = params.pos_range;
		range_sz = range.get_size();
		UNROLL_3X(range_sz_inv[i_] = 1.0/range_sz[i_];)
		clear();
		buildings.reserve(params.num);
		grid.resize(grid_sz*grid_sz); // square
		unsigned num_gen(0);

		for (unsigned i = 0; i < params.num; ++i) {
			building_t b;
			point center;
			
			for (unsigned d = 0; d < 3; ++d) { // x,y,z
				if (d < 2) {center[d] = rgen.rand_uniform(range.d[d][0], range.d[d][1]);} // x,y
				else {center[d] = get_exact_zval(center.x+xlate.x, center.y+xlate.y);} // z
				float const sz(0.5*rgen.rand_uniform(params.sz_range.d [d][0], params.sz_range.d [d][1]));
				b.bcube.d[d][0] = center[d] - ((d == 2) ? 0.0 : sz); // only in XY
				b.bcube.d[d][1] = center[d] + sz;
			} // for d
			if (center.z < water_plane_z) continue; // skip underwater buildings
			++num_gen;

			// check building for overlap with other buildings
			cube_t test_bc(b.bcube);
			test_bc.expand_by(0.1*b.bcube.get_size()); // expand by 10%
			bool overlaps(0);

#if 1 // use grid acceleration
			unsigned ixr[2][2];
			get_grid_range(b.bcube, ixr);

			for (unsigned y = ixr[0][1]; y <= ixr[1][1] && !overlaps; ++y) {
				for (unsigned x = ixr[0][0]; x <= ixr[1][0] && !overlaps; ++x) {
					grid_elem_t const &ge(get_grid_elem(x, y));
					for (auto g = ge.ixs.begin(); g != ge.ixs.end(); ++g) {
						assert(*g < buildings.size());
						if (test_bc.intersects_xy(buildings[*g].bcube)) {overlaps = 1; break;} // Note: only check for XY intersection
					}
				}
			}
#else // brute force iteration
			for (auto j = buildings.begin(); j != buildings.end(); ++j) {
				if (test_bc.intersects_xy(j->bcube)) {overlaps = 1; break;} // Note: only check for XY intersection
			}
#endif
			if (overlaps) continue;
			// FIXME: create acceleration structure: 2D grid in XY plane
			// FIXME: check for overlaps with other buildings
			if (params.flatten_mesh) {flatten_hmap_region(b.bcube);} // flatten the mesh under the bcube to a height of mesh_zval
			for (unsigned d = 0; d < 4; ++d) {b.color[d] = rgen.rand_uniform(params.color_min[d], params.color_max[d]);}
			b.tid = params.tid;
			add_to_grid(b.bcube, buildings.size());
			buildings.push_back(b);
		} // for i
		cout << "Buildings: " << params.num << " / " << num_gen << " / " << buildings.size() << endl;
	}

	void draw(bool shadow_only, vector3d const &xlate) const {
		if (empty()) return;
		fgPushMatrix();
		translate_to(xlate);
		shader_t s;

		if (shadow_only) {
			s.begin_color_only_shader(); // really don't even need colors
		}
		else {
			int const use_bmap(global_building_params.nm_tid >= 0), is_outside(1);
			bool const indir(0), use_smap(0); // FIXME
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, 0, 0, 0, use_smap, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, is_outside);
		}
		for (auto i = buildings.begin(); i != buildings.end(); ++i) {i->draw(s, shadow_only, xlate);}
		s.end_shader();
		fgPopMatrix();
	}
};


building_creator_t building_creator;

void gen_buildings() {building_creator.gen(global_building_params);}
void draw_buildings(bool shadow_only, vector3d const &xlate) {building_creator.draw(shadow_only, xlate);}


