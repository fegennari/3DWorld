// 3D World - OpenGL CS184 Computer Graphics Project - Throwable Sphere Materials
// by Frank Gennari
// 9/5/16

#include "sphere_materials.h"
#include "physics_objects.h"
#include "gameplay.h"
#include "openal_wrap.h"
#include "lightmap.h"
#include "file_utils.h"
#include <fstream>

using namespace std;

unsigned const MAX_SPHERE_MATERIALS = 255;

unsigned spheres_mode(0); // 0=none, 1=spheres, 2=cubes
unsigned max_num_mat_spheres(1);

extern bool spraypaint_mode;
extern int frame_counter;
extern float tfticks, CAMERA_RADIUS, ball_velocity;
extern point sun_pos;
extern colorRGBA sun_color;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern coll_obj_group coll_objects;
extern reflective_cobjs_t reflective_cobjs;
extern vector<light_source_trig> light_sources_d;
extern vector<texture_t> textures;


bool read_texture(char const *const str, unsigned line_num, int &tid, bool is_normal_map, bool invert_y=0);


struct cube_map_lix_t {
	int ixs[6]; // one per cube face, -1 is disabled
	cube_map_lix_t() {for (unsigned i = 0; i < 6; ++i) {ixs[i] = -1;}}
};

class sphere_mat_vect : public vector<sphere_mat_t> {

	unsigned mat_ix;
	typedef map<unsigned, cube_map_lix_t> obj_to_light_map_t;
	obj_to_light_map_t obj_to_light_map;
	vector<unsigned> light_free_list;

public:
	sphere_mat_vect() : mat_ix(0) {}
	unsigned get_ix() const {assert(mat_ix < size()); return mat_ix;}
	sphere_mat_t       &get_cur_mat()       {assert(mat_ix < size()); return operator[](mat_ix);}
	sphere_mat_t const &get_cur_mat() const {assert(mat_ix < size()); return operator[](mat_ix);}
	sphere_mat_t const &get_mat(unsigned ix) const {assert(ix < size()); return operator[](ix);}
	void update_ix(int val) {mat_ix = (mat_ix + size() + val) % size();}

	void remove_light(int light_id) {
		if (light_id < 0) return; // disabled
		assert((unsigned)light_id < light_sources_d.size());
		light_sources_d[light_id].disable();
		light_free_list.push_back(light_id);
	}
	void remove_lights(cube_map_lix_t const &lix) {
		for (unsigned i = 0; i < 6; ++i) {remove_light(lix.ixs[i]);}
	}
	void remove_obj_light(unsigned obj_id) {
		obj_to_light_map_t::iterator it(obj_to_light_map.find(obj_id));
		if (it == obj_to_light_map.end()) return; // not found
		remove_lights(it->second);
		obj_to_light_map.erase(it);
	}
	unsigned alloc_light() {
		if (!light_free_list.empty()) {
			unsigned const light_id(light_free_list.back());
			light_free_list.pop_back();
			return light_id;
		}
		else {
			unsigned const light_id(light_sources_d.size());
			light_sources_d.push_back(light_source_trig());
			return light_id;
		}
	}
	cube_map_lix_t add_obj(unsigned obj_id, bool add_light) {
		remove_obj_light(obj_id);
		cube_map_lix_t ret;
		if (add_light) {
			for (unsigned i = 0; i < 6; ++i) {ret.ixs[i] = alloc_light();}
			obj_to_light_map[obj_id] = ret;
		}
		return ret;
	}
	void sync_light_pos(unsigned obj_id, point const &obj_pos) const {
		obj_to_light_map_t::const_iterator it(obj_to_light_map.find(obj_id));
		if (it == obj_to_light_map.end()) return; // not found
		cube_map_lix_t const &lix(it->second);

		for (unsigned i = 0; i < 6; ++i) {
			unsigned const ix(lix.ixs[i]);
			assert(ix < light_sources_d.size());
			point const old_pos(light_sources_d[ix].get_pos());
			light_sources_d[ix].shift_by(obj_pos - old_pos);
		}
	}
};

sphere_mat_vect sphere_materials;

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
	void read_to_newline() {
		while (1) {
			int const c(in.get());
			if (c == '\n' || c == '\0' || c == EOF) break; // end of file or line
		}
		return;
	}
public:
	material_file_parser_t(string const &fn_) : fn(fn_) {}

	bool read() {
		in.open(fn);
		if (!in.good()) {cerr << "Error: Failed to open sphere materials file '" << fn << "'" << endl; return 0;}
		string key, str;
		sphere_mat_t cur_mat;
		//sphere_materials.clear();

		while (in >> key) {
			if (key[0] == '#') {read_to_newline();} // comment
			else if (key == "add_material") {
				if (!read_mat_value(cur_mat.name, "material name")) return 0;

				if (sphere_materials.size() >= MAX_SPHERE_MATERIALS) {
					cerr << "Error: Too many materials in sphere materials file '" << fn << "': max is " << MAX_SPHERE_MATERIALS << " but saw " << sphere_materials.size() << endl;
					return 0;
				}
				if (cur_mat.tid >= 0 && cur_mat.nm_tid >= 0) { // bind normal map to texture do that material editor can track their relationship
					assert(cur_mat.tid < (int)textures.size() && cur_mat.nm_tid < (int)textures.size());
					textures[cur_mat.tid].maybe_assign_normal_map_tid(cur_mat.nm_tid);
				}
				sphere_materials.push_back(cur_mat);
			}
			else if (key == "shadows") {if (!read_mat_value(cur_mat.shadows, "shadows")) return 0;}
			else if (key == "emissive") {if (!read_mat_value(cur_mat.emissive, "emissive")) return 0;}
			else if (key == "reflective") {if (!read_mat_value(cur_mat.reflective, "reflective")) return 0;}
			else if (key == "destroy_thresh") {if (!read_mat_value(cur_mat.destroy_thresh, "destroy_thresh")) return 0;}
			else if (key == "radius_scale") {if (!read_mat_value(cur_mat.radius_scale, "radius_scale")) return 0;}
			else if (key == "alpha") {if (!read_mat_value(cur_mat.alpha, "alpha")) return 0;}
			else if (key == "metalness") {if (!read_mat_value(cur_mat.metal, "metalness")) return 0;}
			else if (key == "specular_mag") {if (!read_mat_value(cur_mat.spec_mag, "specular_mag")) return 0;}
			else if (key == "specular_exp") {if (!read_mat_value(cur_mat.shine, "specular_exp")) return 0;}
			else if (key == "hardness") {if (!read_mat_value(cur_mat.hardness, "hardness")) return 0;}
			else if (key == "density") {if (!read_mat_value(cur_mat.density, "density")) return 0;}
			else if (key == "refract_ix") {if (!read_mat_value(cur_mat.refract_ix, "refract_ix")) return 0;}
			else if (key == "light_atten") {if (!read_mat_value(cur_mat.light_atten, "light_atten")) return 0;}
			else if (key == "light_radius") {if (!read_mat_value(cur_mat.light_radius, "light_radius")) return 0;}
			else if (key == "diffuse_color") {if (!read_mat_value(cur_mat.diff_c, "diffuse_color")) return 0;}
			else if (key == "specular_color") {if (!read_mat_value(cur_mat.spec_c, "specular_color")) return 0;}
			else if (key == "max_num_spheres") {if (!read_mat_value(max_num_mat_spheres, "max_num_spheres")) return 0;}
			else if (key == "texture") {
				if (!read_mat_value(str, "texture")) return 0;
				if (!read_texture(str.c_str(), 0, cur_mat.tid, 0, 0)) return 0;
			}
			else if (key == "normal_map") {
				if (!read_mat_value(str, "normal_map")) return 0;
				if (!read_texture(str.c_str(), 0, cur_mat.nm_tid, 1, 0)) return 0;
			}
			else {cerr << "Error: Unrecognized keyword in sphere materials file '" << fn << "': " << key << endl; return 0;}
		}
		return 1;
	}
};

bool read_sphere_materials_file(string const &fn) {
	sphere_materials.clear();
	return material_file_parser_t(fn).read();
}

sphere_mat_t &get_cur_sphere_mat() {return sphere_materials.get_cur_mat();}

void show_cur_sphere_mode() {

	if (spheres_mode == 0) {print_text_onscreen("Flashlight", YELLOW, 1.0, TICKS_PER_SECOND, 1); return;}
	string const &str(get_cur_sphere_mat().get_name());
	print_text_onscreen(str, YELLOW, 1.0, TICKS_PER_SECOND, 1); // 1 second
}

void toggle_sphere_mode() {

	if (world_mode != WMODE_GROUND) return;
	if (sphere_materials.empty()) {spheres_mode = 0;} else {spheres_mode = (spheres_mode+1)%3;}
	if (spheres_mode) {spraypaint_mode = 0;}
	show_cur_sphere_mode();
}

void change_sphere_material(int val, bool quiet) {

	if (world_mode != WMODE_GROUND) return;
	sphere_materials.update_ix(val);
	if (quiet) return;
	show_cur_sphere_mode();
	play_switch_weapon_sound();
}


string sphere_mat_t::get_name() const {return name + ((spheres_mode == 2) ? " (Cube)" : " (Sphere)");}

bool throw_sphere(bool mode) {

	static double prev_fticks(0.0);
	if ((double)tfticks - prev_fticks < 20.0) return 0; // 20 ticks = 0.5s fire delay
	prev_fticks = tfticks;

	if (max_num_mat_spheres == 0) return 0;
	int const type(MAT_SPHERE), cid(coll_id[type]);
	if (cid < 0) return 0;
	assert(cid < NUM_TOT_OBJS);
	obj_group &objg(obj_groups[cid]);
	unsigned const mat_ix(sphere_materials.get_ix());
	sphere_mat_t const &mat(sphere_materials.get_mat(mat_ix));
	float const radius(object_types[type].radius*mat.radius_scale), radius_sum(CAMERA_RADIUS + radius);
	int const chosen(objg.choose_object());
	point const fpos(get_camera_pos());
	gen_sound(SOUND_SWING, fpos, 0.5, 1.0);
	objg.create_object_at(chosen, (fpos + cview_dir*radius_sum + plus_z*(0.2*radius_sum)));
	dwobject &obj(objg.get_obj(chosen));
	obj.velocity  = cview_dir*(1.0 + ball_velocity*2.0);
	obj.init_dir  = -cview_dir;
	obj.time      = -1;
	obj.source    = CAMERA_ID;
	if (spheres_mode == 2) {obj.flags |= IS_CUBE_FLAG;}
	assert(mat_ix <= MAX_SPHERE_MATERIALS); // since it's packed into an unsigned char
	obj.direction = (unsigned char)mat_ix;
	bool const has_shadows(mat.light_radius > 0.0 && mat.shadows);
	cube_map_lix_t lix(sphere_materials.add_obj(chosen, has_shadows));

	if (has_shadows) {
		int platform_id(-1); // unused
		float const beamwidth = 0.4; // 0.3 to 0.5 are okay
		float const near_clip(radius);

		for (unsigned ldim = 0; ldim < 3; ++ldim) { // setup 6 light sources, one per cube face
			vector3d dir(zero_vector);

			for (unsigned ldir = 0; ldir < 2; ++ldir) {
				dir[ldim] = (ldir ? 1.0 : -1.0);
				unsigned const ix(lix.ixs[2*ldim + ldir]);
				assert(ix < light_sources_d.size());
				light_source_trig &ls(light_sources_d[ix]);
				ls = light_source_trig(light_source(mat.light_radius, obj.pos, obj.pos, mat.diff_c, 0, dir, beamwidth, 0.0, 1, near_clip), 1, platform_id, 0);
				//ls.bind_to_pos(obj.pos, 1); // dynamic binding
			} // for ldir
		} // for ldim
	}
	return 1;
}

bool is_mat_sphere_a_shadower(dwobject const &obj) {
	sphere_mat_t const &mat(sphere_materials.get_mat(obj.direction));
	//if (mat.light_radius > 0.0) return 0;
	if (mat.alpha < MIN_SHADOW_ALPHA && mat.metal == 0.0) return 0; // transparent glass
	return 1;
}

float get_mat_sphere_density(dwobject const &obj) {return sphere_materials.get_mat(obj.direction).density;}
float get_mat_sphere_rscale (dwobject const &obj) {return sphere_materials.get_mat(obj.direction).radius_scale;}

void sync_mat_sphere_lpos(unsigned id, point const &pos) {sphere_materials.sync_light_pos(id, pos);}

void add_cobj_for_mat_sphere(dwobject &obj, cobj_params const &cp_in) {

	sphere_mat_t const &mat(sphere_materials.get_mat(obj.direction));
	bool const reflective(mat.reflective && enable_all_reflections());
	float const obj_radius(object_types[obj.type].radius*mat.radius_scale); // Note: must match object radius for collision detection to work correctly
	cobj_params cp(cp_in); // deep copy
	cp.draw        = 1; // obj is not drawn
	cp.elastic     = mat.hardness; // elastic is misnamed, really it's hardness
	cp.metalness   = mat.metal;
	cp.is_emissive = mat.emissive;
	cp.color       = colorRGBA(mat.diff_c, mat.alpha);
	cp.spec_color  = mat.spec_c * mat.spec_mag;
	cp.shine       = mat.shine;
	cp.refract_ix  = mat.refract_ix;
	cp.light_atten = mat.light_atten;
	cp.density     = mat.density;
	cp.tscale      = 0.0;
	cp.tid         = mat.tid;
	cp.normal_map  = mat.nm_tid;

	if (obj.flags & IS_CUBE_FLAG) { // cube
		cube_t cube;
		cube.set_from_sphere(obj.pos, obj_radius);
		obj.coll_id = add_coll_cube(cube, cp, -1, 0);
		if (reflective) {add_reflective_cobj(obj.coll_id);}
	}
	else { // sphere
		obj.coll_id = add_coll_sphere(obj.pos, obj_radius, cp, -1, 0, reflective);
	}
	coll_obj &cobj(coll_objects.get_cobj(obj.coll_id));
	cobj.destroy = mat.destroy_thresh;
	
	if (mat.light_radius > 0.0 && !mat.shadows) {
		add_dynamic_light(mat.light_radius, obj.pos, mat.diff_c); // regular point light
	}
	else if (mat.alpha < 0.5 && mat.metal == 0.0 && mat.refract_ix > 1.0) { // glass
		if (is_visible_to_light_cobj(obj.pos, LIGHT_SUN, 0.0, -1, 0)) {
			vector3d const dir((obj.pos - sun_pos).get_norm());
			add_dynamic_light(0.3+obj_radius, obj.pos, sun_color, dir, 0.03); // spotlight to simulate specular focusing; would look better with smoother falloff
		}
	}
	sync_mat_sphere_lpos(cp.cf_index, obj.pos);
}

void remove_mat_sphere(unsigned id) {sphere_materials.remove_obj_light(id);}

