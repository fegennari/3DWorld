// 3D World - OpenGL CS184 Computer Graphics Project - Throwable Sphere Materials
// by Frank Gennari
// 9/5/16

//#include "collision_detect.h"
#include "gameplay.h"
#include "openal_wrap.h"
//#include "shaders.h"
#include <fstream>

using namespace std;

bool spheres_mode(0);
unsigned sphere_material_ix(0);

struct sphere_mat_t {
	bool shadows;
	int destroy_thresh;
	float alpha, metal, spec_mag, shine, hardness;
	colorRGB diff_c, spec_c, emiss_c;
	string name;

	sphere_mat_t() : shadows(0), destroy_thresh(0), alpha(1.0), metal(1.0), spec_mag(0.0), shine(1.0), hardness(0.8), diff_c(WHITE), spec_c(WHITE), emiss_c(BLACK) {}
};

vector<sphere_mat_t> sphere_materials;

class material_file_parser_t {

	string const &fn;
	ifstream in;

	template<typename T> bool read_value(T &val) {return bool(in >> val);}
	bool read_value(colorRGB &val) {return bool(in >> val.R >> val.G >> val.B);}

	template<typename T> bool read_mat_value(T &val, const char *name) {
		if (read_value(val)) return 1;
		cerr << "Error reading " << name << " from sphere materials file '" << fn << "'" << endl;
		return 0;
	}
public:
	material_file_parser_t(string const &fn_) : fn(fn_) {}

	bool read() {
		in.open(fn);
		if (!in.good()) {cerr << "Error: Failed to open sphere materials file '" << fn << "'" << endl; return 0;}
		string key;
		sphere_mat_t cur_mat;

		while (in >> key) {
			if (0) {}
			else if (key == "add_material") {
				if (!read_mat_value(cur_mat.name, "material name")) return 0;
				sphere_materials.push_back(cur_mat);
			}
			else if (key == "shadows") {if (!read_mat_value(cur_mat.shadows, "shadows")) return 0;}
			else if (key == "destroy_thresh") {if (!read_mat_value(cur_mat.destroy_thresh, "destroy_thresh")) return 0;}
			else if (key == "alpha") {if (!read_mat_value(cur_mat.alpha, "alpha")) return 0;}
			else if (key == "metalness") {if (!read_mat_value(cur_mat.metal, "metalness")) return 0;}
			else if (key == "specular_mag") {if (!read_mat_value(cur_mat.spec_mag, "specular_mag")) return 0;}
			else if (key == "specular_exp") {if (!read_mat_value(cur_mat.shine, "specular_exp")) return 0;}
			else if (key == "hardness") {if (!read_mat_value(cur_mat.hardness, "hardness")) return 0;}
			else if (key == "diffuse_color") {if (!read_mat_value(cur_mat.diff_c, "diffuse_color")) return 0;}
			else if (key == "specular_color") {if (!read_mat_value(cur_mat.spec_c, "specular_color")) return 0;}
			else if (key == "emissive_color") {if (!read_mat_value(cur_mat.emiss_c, "emissive_color")) return 0;}
			else {cerr << "Error: Unrecognized keyword in sphere materials file '" << fn << "': " << key << endl; return 0;}
		}
		return 1;
	}
};

bool read_sphere_materials_file(string const &fn) {
	sphere_materials.clear();
	return material_file_parser_t(fn).read();
}

void show_cur_sphere_mode() {

	if (!spheres_mode) {print_text_onscreen("Flashlight", YELLOW, 1.0, TICKS_PER_SECOND, 1); return;}
	assert(sphere_material_ix < sphere_materials.size());
	string const &str(sphere_materials[sphere_material_ix].name);
	print_text_onscreen(str, YELLOW, 1.0, TICKS_PER_SECOND, 1); // 1 second
}

void toggle_sphere_mode() {

	if (world_mode != WMODE_GROUND) return;
	if (sphere_materials.empty()) {spheres_mode = 0;} else {spheres_mode ^= 1;}
	show_cur_sphere_mode();
}

void change_sphere_material(int val) {

	if (world_mode != WMODE_GROUND) return;
	sphere_material_ix = (sphere_material_ix + sphere_materials.size() + val) % sphere_materials.size();
	show_cur_sphere_mode();
	play_switch_weapon_sound();
}

void throw_sphere(bool mode) {

	assert(sphere_material_ix < sphere_materials.size());
	sphere_mat_t const &mat(sphere_materials[sphere_material_ix]);
	point const pos(get_camera_pos());
	// FIXME: WRITE
}
