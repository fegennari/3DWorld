// 3D World - OpenGL CS184 Computer Graphics Project - Throwable Sphere Materials (enable with 'b' + '0')
// by Frank Gennari
// 9/5/16

#include "sphere_materials.h"
#include "physics_objects.h"
#include "openal_wrap.h"
#include "lightmap.h"
#include "file_utils.h"
#include "cube_map_shadow_manager.h"
#include <fstream>

using namespace std;

unsigned const MAX_SPHERE_MATERIALS = 255;
float    const MIN_LIGHT_RADIUS     = 0.01; // very small numbers lead to problems
float    const cube_map_beamwidth   = 0.4; // 0.3 to 0.5 are okay

unsigned spheres_mode(0), num_objs_thrown(0); // 0=none, 1=dynamic spheres, 2=dynamic cubes, 3=static spheres, 4=static cubes
unsigned max_num_mat_spheres(1);
float sphere_mat_fire_delay(0.5); // in seconds

extern bool spraypaint_mode;
extern int frame_counter, display_framerate;
extern float CAMERA_RADIUS, ball_velocity;
extern double tfticks;
extern point sun_pos;
extern colorRGBA sun_color;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern coll_obj_group coll_objects;
extern reflective_cobjs_t reflective_cobjs;
extern vector<light_source> light_sources_a;
extern vector<light_source_trig> light_sources_d;
extern vector<texture_t> textures;


void cube_map_lix_t::add_cube_face_lights(point const &pos, float radius, colorRGBA const &color, float near_clip, bool outdoor_shadows) {

	for (unsigned ldim = 0; ldim < 3; ++ldim) { // setup 6 light sources, one per cube face
		vector3d dir(zero_vector);
		for (unsigned ldir = 0; ldir < 2; ++ldir) {
			dir[ldim] = (ldir ? 1.0 : -1.0);
			unsigned const ix(ixs[2*ldim + ldir]);
			assert(ix < light_sources_d.size());
			light_source_trig &ls(light_sources_d[ix]);
			ls = light_source_trig(light_source(radius, pos, pos, color, 0, dir, cube_map_beamwidth, 0.0, 1, near_clip), 1, -1, 0, sensor_t(), outdoor_shadows);
			//ls.bind_to_pos(obj.pos, 1); // dynamic binding
		} // for ldir
	} // for ldim
}

void cube_map_shadow_manager::remove_light(int light_id) {
	if (light_id < 0) return; // disabled
	assert((unsigned)light_id < light_sources_d.size());
	light_sources_d[light_id].disable();
	light_free_list.push_back(light_id);
}
void cube_map_shadow_manager::remove_lights(cube_map_lix_t const &lix) {
	for (unsigned i = 0; i < 6; ++i) {remove_light(lix.ixs[i]);}
}
void cube_map_shadow_manager::remove_obj_light(unsigned obj_id) {
	obj_to_light_map_t::iterator it(obj_to_light_map.find(obj_id));
	if (it == obj_to_light_map.end()) return; // not found
	remove_lights(it->second);
	obj_to_light_map.erase(it);
}
unsigned cube_map_shadow_manager::alloc_light() {
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
cube_map_lix_t cube_map_shadow_manager::add_obj(unsigned obj_id, bool add_light) {
	remove_obj_light(obj_id);
	cube_map_lix_t ret;
	
	if (add_light) {
		for (unsigned i = 0; i < 6; ++i) {ret.ixs[i] = alloc_light();}
		obj_to_light_map[obj_id] = ret;
	}
	return ret;
}
void cube_map_shadow_manager::sync_light_pos(unsigned obj_id, point const &obj_pos) const {
	obj_to_light_map_t::const_iterator it(obj_to_light_map.find(obj_id));
	if (it == obj_to_light_map.end()) return; // not found
	cube_map_lix_t const &lix(it->second);

	for (unsigned i = 0; i < 6; ++i) {
		unsigned const ix(lix.ixs[i]);
		assert(ix < light_sources_d.size());
		light_sources_d[ix].move_to(obj_pos);
	}
}


bool read_texture(char const *const str, unsigned line_num, int &tid, bool is_normal_map, bool invert_y=0);


class sphere_mat_vect : public vector<sphere_mat_t>, public cube_map_shadow_manager {

	unsigned mat_ix;

public:
	sphere_mat_vect() : mat_ix(0) {}
	unsigned get_ix() const {assert(mat_ix < size()); return mat_ix;}
	sphere_mat_t       &get_cur_mat()       {assert(mat_ix < size()); return operator[](mat_ix);}
	sphere_mat_t const &get_cur_mat() const {assert(mat_ix < size()); return operator[](mat_ix);}
	sphere_mat_t const &get_mat(unsigned ix) const {assert(ix < size()); return operator[](ix);}
	void update_ix(int val) {mat_ix = (mat_ix + size() + val) % size();}

	static string texture_str(int tid) {
		if (tid < 0) return "none";
		assert((unsigned)tid < textures.size());
		return textures[tid].name;
	}
	bool write_to_file(string const &fn) const {
		ofstream out(fn);
		if (!out.good()) {cerr << "Error: Failed to open sphere materials file '" << fn << "' for writing" << endl; return 0;}
		out << "max_num_spheres " << max_num_mat_spheres << endl;
		out << "fire_delay " << sphere_mat_fire_delay << endl;

		for (const_iterator i = begin(); i != end(); ++i) {
			out << endl;
			out << "shadows " << i->shadows << endl;
			out << "emissive " << i->emissive << endl;
			out << "reflective " << i->reflective << endl;
			out << "destroyable " << i->destroyable << endl;
			out << "radius_scale " << i->radius_scale << endl;
			out << "alpha " << i->alpha << endl;
			out << "metalness " << i->metal << endl;
			out << "specular_mag " << i->spec_mag << endl;
			out << "specular_exp " << i->shine << endl;
			out << "hardness " << i->hardness << endl;
			out << "density " << i->density << endl;
			out << "refract_ix " << i->refract_ix << endl;
			out << "light_atten " << i->light_atten << endl;
			out << "light_radius " << i->light_radius << endl;
			out << "diffuse_color " << i->diff_c.raw_str() << endl;
			out << "specular_color " << i->spec_c.raw_str() << endl;
			out << "texture " << texture_str(i->tid) << endl;
			out << "normal_map " << texture_str(i->nm_tid) << endl;
			out << "add_material " << i->name << endl;
		}
		return 1;
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
		cerr << "Error reading " << name << " from sphere materials file '" << fn << "' for reading" << endl;
		return 0;
	}
	void read_to_newline() {
		while (1) {
			int const c(in.get());
			if (c == '\n' || is_EOF(c)) break; // end of file or line
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
				if (cur_mat.tid >= 0 && cur_mat.nm_tid >= 0) { // bind normal map to texture so that material editor can track their relationship
					assert(cur_mat.tid < (int)textures.size() && cur_mat.nm_tid < (int)textures.size());
					textures[cur_mat.tid].maybe_assign_normal_map_tid(cur_mat.nm_tid);
				}
				sphere_materials.push_back(cur_mat);
			}
			else if (key == "fire_delay") {if (!read_mat_value(sphere_mat_fire_delay, "fire_delay")) return 0;}
			else if (key == "shadows") {if (!read_mat_value(cur_mat.shadows, "shadows")) return 0;}
			else if (key == "emissive") {if (!read_mat_value(cur_mat.emissive, "emissive")) return 0;}
			else if (key == "reflective") {if (!read_mat_value(cur_mat.reflective, "reflective")) return 0;}
			else if (key == "destroyable") {if (!read_mat_value(cur_mat.destroyable, "destroyable")) return 0;} // only {0=none, 1=shatterable, 2=explodeable}
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
	if (fn.empty()) return 1;
	sphere_materials.clear();
	return material_file_parser_t(fn).read();
}
bool write_sphere_materials_file(string const &fn) {
	if (fn.empty()) return 1;
	return sphere_materials.write_to_file(fn);
}

sphere_mat_t &get_cur_sphere_mat() {return sphere_materials.get_cur_mat();}

void show_cur_sphere_mode() {

	if (spheres_mode == 0) {print_text_onscreen("Flashlight", YELLOW, 1.0, TICKS_PER_SECOND, 1); return;}
	string const &str(get_cur_sphere_mat().get_name());
	print_text_onscreen(str, YELLOW, 1.0, TICKS_PER_SECOND, 1); // 1 second
}

void toggle_sphere_mode() {

	if (world_mode != WMODE_GROUND) return;
	if (sphere_materials.empty()) {spheres_mode = 0;} else {spheres_mode = (spheres_mode+1)%5;}
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


string const mode_strs[5] = {"None", "Dynamic Sphere", "Dynamic Cube", "Static Sphere", "Static Cube"};

string sphere_mat_t::get_name() const {return name + " (" + mode_strs[spheres_mode] + ")";}


void set_cobj_params_from_material(cobj_params &cp, sphere_mat_t const &mat) {

	cp.draw        = 1; // obj is drawn
	cp.elastic     = mat.hardness; // elastic is misnamed, really it's hardness
	cp.metalness   = mat.metal;
	cp.is_emissive = mat.emissive;
	cp.color       = colorRGBA(mat.diff_c, mat.alpha);
	cp.spec_color  = mat.spec_c * mat.spec_mag;
	cp.shine       = mat.shine;
	cp.refract_ix  = mat.refract_ix;
	cp.light_atten = mat.light_atten;
	cp.density     = mat.density;
	cp.tscale      = ((mat.tid >= 0 || mat.alpha == 1.0) ? 0.0 : 1.0); // use tscale for untextured transparent materials for back-to-front draw order in coll_obj::draw_cobj()
	cp.tid         = mat.tid;
	cp.normal_map  = mat.nm_tid;
}

int add_cobj_with_material(cobj_params const &cp, sphere_mat_t const &mat, point const &pos, float base_radius, bool is_cube, bool is_static) {

	bool const reflective(mat.reflective && (is_static || enable_all_reflections()));
	float const obj_radius(base_radius*mat.radius_scale);
	int coll_id(-1);

	if (is_cube) { // cube
		cube_t cube;
		cube.set_from_sphere(pos, obj_radius);
		coll_id = add_coll_cube(cube, cp, -1, 0);
	}
	else { // sphere
		coll_id = add_coll_sphere(pos, obj_radius, cp, -1, 0);
		//coll_id = add_coll_torus(pos, /*signed_rand_vector_norm()*/plus_z, obj_radius, 0.4*obj_radius, cp, -1, 0);
		//coll_id = add_coll_cylinder(pos, (pos + signed_rand_vector(obj_radius*rand_uniform(0.5, 2.0))), obj_radius, obj_radius*rand_uniform(0.5, 1.0), cp, -1, 0);
		//coll_id = add_coll_cylinder(pos, (pos + vector3d(0.0, 0.0, 2.0*obj_radius)), obj_radius, obj_radius, cp, -1, 0);
	}
	if (reflective) {add_reflective_cobj(coll_id);}
	if (is_static) {coll_objects.get_cobj(coll_id).destroy = 2*min(2, max(0, mat.destroyable));} // only {0=none, 1=shatterable, 2=explodeable}
	return coll_id;
}


void add_static_material_object(cobj_params &cp, sphere_mat_t const &mat, point const &pos, float base_radius, bool is_cube) {

	bool const has_shadows(mat.light_radius > MIN_LIGHT_RADIUS && mat.shadows);
	cp.flags |= COBJ_MOVABLE;
	set_cobj_params_from_material(cp, mat);
	int const coll_id(add_cobj_with_material(cp, mat, pos, base_radius, is_cube, 1));
	assert(coll_id >= 0);
	coll_obj &cobj(coll_objects.get_cobj(coll_id));
	cobj.fixed = 1; // make it static
	register_moving_cobj(coll_id); // let it fall
	register_movable_cobj_shadow(coll_id);
	vector<light_source> lss;
	float const sphere_radius(base_radius*mat.radius_scale);

	if (has_shadows) { // special case for shadowed point light with cube map
		float const near_clip(1.01*sphere_radius);

		for (unsigned ldim = 0; ldim < 3; ++ldim) {
			vector3d dir(zero_vector);
			for (unsigned ldir = 0; ldir < 2; ++ldir) {
				dir[ldim] = (ldir ? 1.0 : -1.0);
				lss.push_back(light_source(mat.light_radius, pos, pos, mat.diff_c, 0, dir, cube_map_beamwidth, sphere_radius, 1, near_clip)); // add lights for each cube face
			}
		}
	}
	else if (mat.light_radius > MIN_LIGHT_RADIUS) {
		lss.push_back(light_source(mat.light_radius, pos, pos, mat.diff_c, 0, zero_vector, 1.0, sphere_radius));
	}
	for (auto ls = lss.begin(); ls != lss.end(); ++ls) {
		light_sources_d.push_back(light_source_trig(*ls, has_shadows, -1, 0));
		light_sources_d.back().bind_to_pos(pos, 0, coll_id);
	}
}


bool throw_sphere(bool mode) {

	static double prev_fticks(0.0);
	if ((tfticks - prev_fticks) < sphere_mat_fire_delay*double(TICKS_PER_SECOND)) return 0; // 20 ticks = 0.5s fire delay
	prev_fticks = tfticks;

	if (max_num_mat_spheres == 0) return 0;
	int const type(MAT_SPHERE), cid(coll_id[type]);
	if (cid < 0) return 0;
	unsigned const mat_ix(sphere_materials.get_ix());
	sphere_mat_t const &mat(sphere_materials.get_mat(mat_ix));
	float const base_radius(object_types[type].radius), radius(base_radius*mat.radius_scale), radius_sum(CAMERA_RADIUS + radius), near_clip(1.01*radius);
	bool const is_cube(spheres_mode == 2 || spheres_mode == 4);
	point const fpos(get_camera_pos() + cview_dir*radius_sum*(is_cube ? SQRT2 : 1.0) + plus_z*(0.2*radius_sum));
	++num_objs_thrown;
	gen_sound(SOUND_SWING, fpos, 0.5, 1.0);

	if (spheres_mode == 3 || spheres_mode == 4) { // static objects
		cobj_params cp(0.0, BLACK, 1, 0); // elastic, color, draw, is_dynamic
		add_static_material_object(cp, mat, fpos, base_radius, is_cube);
		if (display_framerate) {print_debug_text(std::to_string(num_objs_thrown), -100);} // low priority
		return 1;
	}
	assert(cid < NUM_TOT_OBJS);
	obj_group &objg(obj_groups[cid]);
	int const chosen(objg.choose_object());
	objg.create_object_at(chosen, fpos);
	dwobject &obj(objg.get_obj(chosen));
	obj.velocity = cview_dir*(1.0 + ball_velocity*2.0);
	obj.init_dir = -cview_dir;
	obj.time     = -1;
	obj.source   = CAMERA_ID;
	if (is_cube) {obj.flags |= IS_CUBE_FLAG;}
	assert(mat_ix <= MAX_SPHERE_MATERIALS); // since it's packed into an unsigned char
	obj.direction = (unsigned char)mat_ix;
	bool const has_shadows(mat.light_radius > MIN_LIGHT_RADIUS && mat.shadows);
	cube_map_lix_t lix(sphere_materials.add_obj(chosen, has_shadows));
	if (has_shadows) {lix.add_cube_face_lights(obj.pos, mat.light_radius, mat.diff_c, near_clip);}
	return 1;
}

bool is_mat_sphere_a_shadower(dwobject const &obj) {
	sphere_mat_t const &mat(sphere_materials.get_mat(obj.direction));
	//if (mat.light_radius > MIN_LIGHT_RADIUS) return 0;
	if (mat.alpha < MIN_SHADOW_ALPHA && mat.metal == 0.0) return 0; // transparent glass
	return 1;
}

float get_mat_sphere_density(dwobject const &obj) {return sphere_materials.get_mat(obj.direction).density;}
float get_mat_sphere_rscale (dwobject const &obj) {return sphere_materials.get_mat(obj.direction).radius_scale;}
void sync_mat_sphere_lpos(unsigned id, point const &pos) {sphere_materials.sync_light_pos(id, pos);}

void add_cobj_for_mat_sphere(dwobject &obj, cobj_params const &cp_in) {

	sphere_mat_t const &mat(sphere_materials.get_mat(obj.direction));
	float const base_radius(object_types[obj.type].radius); // Note: must match object radius for collision detection to work correctly
	bool const is_cube((obj.flags & IS_CUBE_FLAG) != 0);
	cobj_params cp(cp_in); // deep copy
	set_cobj_params_from_material(cp, mat);
	obj.coll_id = add_cobj_with_material(cp, mat, obj.pos, base_radius, is_cube, 0);
	coll_obj &cobj(coll_objects.get_cobj(obj.coll_id));
	
	if (mat.light_radius > MIN_LIGHT_RADIUS && !mat.shadows) {
		add_dynamic_light(mat.light_radius, obj.pos, mat.diff_c); // regular point light
	}
	else if (mat.alpha < 0.5 && mat.metal == 0.0 && mat.refract_ix > 1.0 && !is_cube) { // glass
		if (is_visible_to_light_cobj(obj.pos, LIGHT_SUN, 0.0, -1, 0)) {
			vector3d const dir((obj.pos - sun_pos).get_norm());
			add_dynamic_light(0.3+cobj.radius, obj.pos, sun_color, dir, 0.03); // spotlight to simulate specular focusing; would look better with smoother falloff
		}
	}
	sync_mat_sphere_lpos(cp.cf_index, obj.pos);
}

void remove_mat_sphere(unsigned id) {sphere_materials.remove_obj_light(id);}


struct gen_sphere_params_t {

	bool enable_reflect, enable_transparent, enable_light_atten, enable_shadows;
	float metal_prob, emissive_prob, metal_white_prob, emiss_white_prob, max_light_atten, max_light_radius;
	int rand_seed;

	gen_sphere_params_t() : enable_reflect(1), enable_transparent(1), enable_light_atten(1), enable_shadows(1),
		metal_prob(0.2), emissive_prob(0.25), metal_white_prob(0.5), emiss_white_prob(0.5), max_light_atten(20.0), max_light_radius(10.0), rand_seed(0) {}
	static bool read_error(string const &str) {cout << "Error reading sphere_gen config option " << str << "." << endl; return 0;}

	bool read_option(FILE *fp) {
		char strc[MAX_CHARS] = {0};
		if (!read_str(fp, strc)) return 0;
		string const str(strc);

		if (str == "enable_reflect") {
			if (!read_bool(fp, enable_reflect)) {return read_error(str);}
		}
		else if (str == "enable_transparent") {
			if (!read_bool(fp, enable_transparent)) {return read_error(str);}
		}
		else if (str == "enable_light_atten") {
			if (!read_bool(fp, enable_light_atten)) {return read_error(str);}
		}
		else if (str == "enable_shadows") {
			if (!read_bool(fp, enable_shadows)) {return read_error(str);}
		}
		else if (str == "metal_prob") {
			if (!read_float(fp, metal_prob) || metal_prob < 0.0 || metal_prob > 1.0) {return read_error(str);}
		}
		else if (str == "emissive_prob") {
			if (!read_float(fp, emissive_prob) || emissive_prob < 0.0 || emissive_prob > 1.0) {return read_error(str);}
		}
		else if (str == "metal_white_prob") {
			if (!read_float(fp, metal_white_prob) || metal_white_prob < 0.0 || metal_white_prob > 1.0) {return read_error(str);}
		}
		else if (str == "emiss_white_prob") {
			if (!read_float(fp, emiss_white_prob) || emiss_white_prob < 0.0 || emiss_white_prob > 1.0) {return read_error(str);}
		}
		else if (str == "max_light_atten") {
			if (!read_float(fp, max_light_atten) || max_light_atten < 0.0) {return read_error(str);}
		}
		else if (str == "max_light_radius") {
			if (!read_float(fp, max_light_radius) || max_light_radius < 0.0) {return read_error(str);}
		}
		else if (str == "rand_seed") {
			if (!read_int(fp, rand_seed)) {return read_error(str);}
		}
		else {
			cout << "Unrecognized sphere_gen keyword in input file: " << str << endl;
			return 0;
		}
		return 1;
	}
};

gen_sphere_params_t gen_sphere_params;

bool parse_sphere_gen_option(FILE *fp) {return gen_sphere_params.read_option(fp);}

// Note: update reflections with 'i'
void gen_rand_spheres(unsigned num, point const &center, float place_radius, float min_radius, float max_radius) {

	timer_t timer("Gen Rand Spheres");
	rand_gen_t rgen;
	vector<sphere_t> spheres;
	gen_sphere_params_t const &sp(gen_sphere_params);
	if (sp.rand_seed != 0) {rgen.set_state(sp.rand_seed, 123);}

	for (unsigned n = 0; n < num; ++n) {
		float const radius(rgen.rand_uniform(min_radius, max_radius));
		vector3d v;
		point pos;

		for (unsigned N = 0; N < 1000; ++N) { // make 100 placement attempts
			while (1) {
				v = rgen.signed_rand_vector_xy();
				if (v.mag_sq() < 1.0) break;
			}
			pos = center + place_radius*v + vector3d(0.0, 0.0, radius);
			bool overlap(0);

			for (auto i = spheres.begin(); i != spheres.end(); ++i) { // check for overlap with previously placed spheres
				if (dist_less_than(pos, i->pos, (radius + i->radius))) {overlap = 1; break;}
			}
			if (!overlap) break;
		}
		cobj_params cp(0.0, BLACK, 1, 0); // elastic, color, draw, is_dynamic
		spheres.emplace_back(pos, radius);
		sphere_mat_t mat;
		bool const is_metal(sp.enable_reflect && rgen.rand_float() < sp.metal_prob);
		mat.metal      = (is_metal ? 1.0 : 0.0);
		mat.spec_mag   = (is_metal ? 1.0 : CLIP_TO_01(rgen.rand_uniform(-0.5, 1.2)));
		mat.shine      = rgen.rand_uniform(1.0, 8.0)*rgen.rand_uniform(1.0, 8.0); // or roughness for reflective objects/metals; 1.0-64.0
		mat.reflective = (sp.enable_reflect && mat.spec_mag > 0.75); // metal is always reflective
		mat.emissive   = (!mat.reflective && rgen.rand_float() < sp.emissive_prob);
		if (!mat.emissive && !is_metal && sp.enable_transparent) {mat.alpha = CLIP_TO_01(rgen.rand_uniform((mat.reflective ? -2.0 : 0.25), 2.0));}
		mat.shadows    = (sp.enable_shadows && mat.alpha > 0.5);
		mat.density    = (is_metal ? 2.0 : 1.0)*rgen.rand_uniform(0.5, 4.0);
		if (sp.max_light_atten > 0.0 && mat.alpha < 0.5) {mat.light_atten = max(rgen.rand_uniform(-sp.max_light_atten, sp.max_light_atten), 0.0f);}
		mat.refract_ix   = rgen.rand_uniform(1.0, 1.5)*rgen.rand_uniform(1.0, 1.5)*rgen.rand_uniform(1.0, 1.5); // 1.0-6.25
		if (sp.max_light_radius > 0.0 && mat.emissive) {mat.light_radius = rgen.rand_uniform(0.5*sp.max_light_radius, 1.0*sp.max_light_radius)*radius;}
		colorRGBA color; // mat.alpha?
		if (is_metal && rgen.rand_float() < sp.metal_white_prob) {color = WHITE;} // make some metals white
		else if (mat.light_radius > 0.0 && rgen.rand_float() < sp.emiss_white_prob) {color = WHITE;} // make some lights white
		else {
			for (unsigned i = 0; i < 3; ++i) {color[i] = CLIP_TO_01(rgen.rand_uniform(-0.25, 1.5));} // saturate in some cases
		}
		if (is_metal) {
			mat.diff_c = BLACK;
			mat.spec_c = color;
		}
		else {
			mat.diff_c = color;
			mat.spec_c = WHITE;
		}
		set_cobj_params_from_material(cp, mat);
		add_static_material_object(cp, mat, pos, radius, 0); // is_cube=0
		
		if (mat.light_radius > MIN_LIGHT_RADIUS) { // needed for indir lighting
			light_sources_a.push_back(light_source(mat.light_radius, pos, pos, mat.diff_c, 0, zero_vector, 1.0, 1.01*radius)); // slightly larger than radius
		}
	} // for n
}

