// 3D World - Building Generation
// by Frank Gennari
// 5/22/17

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "file_utils.h"

using std::string;


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

	void draw(bool shadow_only, vector3d const &xlate=zero_vector) const { // FIXME: more efficient, use batched verts or VBO
		if (!camera_pdu.sphere_and_cube_visible_test((bcube.get_cube_center() + xlate), bcube.get_bsphere_radius(), (bcube + xlate))) return; // FIXME: cache center and radius?
		if (!shadow_only) {select_texture(tid); select_multitex(nm_tid, 5);}
		draw_simple_cube(bcube, (!shadow_only && tid >= 0));
	}
};

class building_creator_t {

	rand_gen_t rgen;
	vector<building_t> buildings;
public:
	bool empty() const {return buildings.empty();}
	void gen(building_params_t const &params) {
		timer_t timer("Gen Buildings");
		buildings.resize(params.num);

		for (unsigned i = 0; i < params.num; ++i) {
			building_t &b(buildings[i]);
			point center;
			
			for (unsigned d = 0; d < 3; ++d) { // x,y,z
				if (d < 2) { // x,y
					center[d] = rgen.rand_uniform(params.pos_range.d[d][0], params.pos_range.d[d][1]);
				}
				else { // z
					center[d] = get_exact_zval(center.x, center.y);
				}
				float const sz(0.5*rgen.rand_uniform(params.sz_range.d [d][0], params.sz_range.d [d][1]));
				b.bcube.d[d][0] = center[d] - sz;
				b.bcube.d[d][1] = center[d] + sz;
			} // for d
			// FIXME: check for overlaps with other buildings
			// FIXME: create acceleration structure: 2D grid in XY plane
			if (params.flatten_mesh) {flatten_hmap_region(b.bcube);} // flatten the mesh under the bcube to a height of mesh_zval
			for (unsigned d = 0; d < 4; ++d) {b.color[d] = rgen.rand_uniform(params.color_min[d], params.color_max[d]);}
			b.tid = params.tid;
		} // for i
	}

	void draw(bool shadow_only, vector3d const &xlate) const {
		if (empty()) return;
		shader_t s;

		if (shadow_only) {
			s.begin_color_only_shader(); // really don't even need colors
		}
		else {
			int const use_bmap(global_building_params.nm_tid >= 0), is_outside(1);
			bool const indir(0), use_smap(0); // FIXME
			setup_smoke_shaders(s, 0.0, 0, 0, indir, 1, 0, 0, 0, use_smap, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, is_outside);
		}
		for (auto i = buildings.begin(); i != buildings.end(); ++i) {i->draw(shadow_only, xlate);}
		s.end_shader();
	}
};


building_creator_t building_creator;

void gen_buildings() {building_creator.gen(global_building_params);}
void draw_buildings(bool shadow_only, vector3d const &xlate) {building_creator.draw(shadow_only, xlate);}


