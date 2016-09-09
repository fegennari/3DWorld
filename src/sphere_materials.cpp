// 3D World - OpenGL CS184 Computer Graphics Project - Throwable Sphere Materials
// by Frank Gennari
// 9/5/16

#include "physics_objects.h"
#include "gameplay.h"
#include "openal_wrap.h"
#include <fstream>

using namespace std;

unsigned const MAX_SPHERE_MATERIALS = 255;

bool spheres_mode(0);
unsigned sphere_material_ix(0), max_num_mat_spheres(1);

extern int frame_counter;
extern float tfticks, CAMERA_RADIUS, ball_velocity;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern coll_obj_group coll_objects;
extern reflective_cobjs_t reflective_cobjs;


struct sphere_mat_t {
	bool light, shadows, emissive, reflective;
	int destroy_thresh;
	float alpha, metal, spec_mag, shine, hardness, density, refract_ix, light_radius;
	colorRGB diff_c, spec_c;
	string name;

	sphere_mat_t() : light(0), shadows(0), emissive(0), reflective(0), destroy_thresh(0), alpha(1.0), metal(1.0),
		spec_mag(0.0), shine(1.0), hardness(0.8), density(1.0), refract_ix(1.0), light_radius(0.0), diff_c(WHITE), spec_c(WHITE) {}
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
		//sphere_materials.clear();

		while (in >> key) {
			if (0) {}
			else if (key == "add_material") {
				if (!read_mat_value(cur_mat.name, "material name")) return 0;

				if (sphere_materials.size() >= MAX_SPHERE_MATERIALS) {
					cerr << "Error: Too many materials in sphere materials file '" << fn << "': max is " << MAX_SPHERE_MATERIALS << " but saw " << sphere_materials.size() << endl;
					return 0;
				}
				sphere_materials.push_back(cur_mat);
			}
			else if (key == "light") {if (!read_mat_value(cur_mat.light, "light")) return 0;}
			else if (key == "shadows") {if (!read_mat_value(cur_mat.shadows, "shadows")) return 0;}
			else if (key == "emissive") {if (!read_mat_value(cur_mat.emissive, "emissive")) return 0;}
			else if (key == "reflective") {if (!read_mat_value(cur_mat.reflective, "reflective")) return 0;}
			else if (key == "destroy_thresh") {if (!read_mat_value(cur_mat.destroy_thresh, "destroy_thresh")) return 0;}
			else if (key == "alpha") {if (!read_mat_value(cur_mat.alpha, "alpha")) return 0;}
			else if (key == "metalness") {if (!read_mat_value(cur_mat.metal, "metalness")) return 0;}
			else if (key == "specular_mag") {if (!read_mat_value(cur_mat.spec_mag, "specular_mag")) return 0;}
			else if (key == "specular_exp") {if (!read_mat_value(cur_mat.shine, "specular_exp")) return 0;}
			else if (key == "hardness") {if (!read_mat_value(cur_mat.hardness, "hardness")) return 0;}
			else if (key == "density") {if (!read_mat_value(cur_mat.density, "density")) return 0;}
			else if (key == "refract_ix") {if (!read_mat_value(cur_mat.refract_ix, "refract_ix")) return 0;}
			else if (key == "light_radius") {if (!read_mat_value(cur_mat.light_radius, "light_radius")) return 0;}
			else if (key == "diffuse_color") {if (!read_mat_value(cur_mat.diff_c, "diffuse_color")) return 0;}
			else if (key == "specular_color") {if (!read_mat_value(cur_mat.spec_c, "specular_color")) return 0;}
			else if (key == "max_num_spheres") {if (!read_mat_value(max_num_mat_spheres, "max_num_spheres")) return 0;}
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

bool throw_sphere(bool mode) {

	static double prev_fticks(0.0);
	if ((double)tfticks - prev_fticks < 20.0) return 0; // 20 ticks = 0.5s fire delay
	prev_fticks = tfticks;

	if (max_num_mat_spheres == 0) return 0;
	int const type(MAT_SPHERE), cid(coll_id[type]);
	if (cid < 0) return 0;
	assert(cid < NUM_TOT_OBJS);
	obj_group &objg(obj_groups[cid]);
	float const radius_sum(CAMERA_RADIUS + object_types[type].radius);
	int const chosen(objg.choose_object());
	point const fpos(get_camera_pos());
	gen_sound(SOUND_SWING, fpos, 0.5, 1.0);
	objg.create_object_at(chosen, (fpos + cview_dir*radius_sum + plus_z*(0.2*radius_sum)));
	dwobject &obj(objg.get_obj(chosen));
	obj.velocity  = cview_dir*(1.0 + ball_velocity*2.0);
	obj.init_dir  = -cview_dir;
	obj.time      = -1;
	obj.source    = CAMERA_ID;

	assert(sphere_material_ix <= MAX_SPHERE_MATERIALS); // since it's packed into an unsigned char
	assert(sphere_material_ix < sphere_materials.size());
	obj.direction = (unsigned char)sphere_material_ix;
	return 1;
}

void add_cobj_for_mat_sphere(dwobject &obj, cobj_params const &cp_in) {

	unsigned const sphere_material_ix(obj.direction);
	assert(sphere_material_ix < sphere_materials.size());
	sphere_mat_t const &mat(sphere_materials[sphere_material_ix]);
	bool const reflective(mat.reflective && enable_all_reflections());
	float const obj_radius(object_types[obj.type].radius); // Note: must match object radius for collision detection to work correctly
	float const light_radius((mat.light_radius == 0.0) ? 20.0*obj_radius : mat.light_radius); // use default radius if material radius is zero
	cobj_params cp(cp_in); // deep copy
	cp.draw        = 1; // obj is not drawn
	cp.elastic     = mat.hardness; // elastic is misnamed, really it's hardness
	cp.metalness   = mat.metal;
	cp.is_emissive = mat.emissive;
	cp.color       = colorRGBA(mat.diff_c, mat.alpha);
	cp.spec_color  = mat.spec_c * mat.spec_mag;
	cp.shine       = mat.shine;
	cp.refract_ix  = mat.refract_ix;
	cp.density     = mat.density;
	cp.tscale      = 0.0;
	cp.tid         = -1;
	obj.coll_id    = add_coll_sphere(obj.pos, obj_radius, cp, -1, 0, reflective);
	coll_obj &cobj(coll_objects.get_cobj(obj.coll_id));
	cobj.destroy   = mat.destroy_thresh;
	
	if (mat.light) {
		add_dynamic_light(light_radius, obj.pos, mat.diff_c);
		if (mat.shadows) {} // FIXME
	}
}
