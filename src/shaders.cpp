// 3D World - Vertex and Fragment GLSL Shader Framework
// by Frank Gennari
// 11/1/10
#include "shaders.h"
#include "mesh.h" // for scene bounds
#include "gl_ext_arb.h"
#include "transform_obj.h"
#include <fstream>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>

using namespace std;

bool const PRINT_SHADER = 0; // print shaders loaded from files
bool const DEBUG_SHADER = 0; // print final generated shaders
bool const PRINT_LOG    = 0;
bool const GEN_FINAL_SHADER_FILES = 0;

string const shaders_dir = "shaders";
string const shader_prefix_shared_file = "common_header_shared.part*";
string const shader_name_table  [NUM_SHADER_TYPES] = {"vert", "frag", "geom", "tesc", "tese", "comp"};
string const shader_prefix_files[NUM_SHADER_TYPES] = {"common_header", "common_header", "common_header", "", "", ""}; // always included
// future work (Nvidia Turing only): GL_MESH_SHADER_NV, GL_TASK_SHADER_NV
int const shader_type_table[NUM_SHADER_TYPES] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_GEOMETRY_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER, GL_COMPUTE_SHADER};

shader_t *cur_shader(NULL);

extern bool fog_enabled, mvm_changed, init_core_context;
extern int is_cloudy, display_mode, fullscreen;
extern unsigned enabled_lights;
extern float cur_fog_end;
extern colorRGBA cur_fog_color;
extern gl_light_params_t gl_light_params[MAX_SHADER_LIGHTS];


void set_one_texture(shader_t &s, unsigned tid, unsigned tu_id, const char *const name) {
	bind_texture_tu(tid, tu_id);
	s.add_uniform_int(name, tu_id);
}

string property_map_t::empty_str;


// *** uniform variables setup ***


char const *append_ix(string &s, unsigned i, bool as_array) {

	assert(i <= 9);
	if (as_array) {s.push_back('[');}
	s.push_back('0'+i);
	if (as_array) {s.push_back(']');}
	return s.c_str();
}


int shader_t::get_uniform_loc(char const *const name) const {

	assert(program && name);
	int const loc(glGetUniformLocation(program, name));
	//cout << "name: " << name << ", loc: " << loc << endl;
	//assert(loc >= 0); // Note: if variable is unused, loc will be -1
	return loc;
}

// this block is all static member functions
bool shader_t::set_uniform_float_array(int loc, float const *const val, unsigned num) {
	if (loc >= 0) {glUniform1fv(loc, num, val); return 1;} else {return 0;}
}
bool shader_t::set_uniform_float(int loc, float val) {
	if (loc >= 0) {glUniform1f(loc, val); return 1;} else {return 0;}
}
bool shader_t::set_uniform_int(int loc, int val) {
	if (loc >= 0) {glUniform1i(loc, val); return 1;} else {return 0;}
}
bool shader_t::set_uniform_uint(int loc, unsigned val) {
	if (loc >= 0) {glUniform1ui(loc, val); return 1;} else {return 0;}
}
bool shader_t::set_uniform_handle(int loc, GLuint64 val) {
	if (loc >= 0) {glUniformHandleui64ARB(loc, val); return 1;} else {return 0;}
}
bool shader_t::set_uniform_vector2d(int loc, vector2d const &val) {
	if (loc >= 0) {glUniform2fv(loc, 1, &val.x); return 1;} else {return 0;}
}
bool shader_t::set_uniform_vector3d(int loc, vector3d const &val) { // same as colorRGB
	if (loc >= 0) {glUniform3fv(loc, 1, &val.x); return 1;} else {return 0;}
}
bool shader_t::set_uniform_vector4d(int loc, vector4d const &val) { // same as colorRGBA
	if (loc >= 0) {glUniform4fv(loc, 1, &val.x); return 1;} else {return 0;}
}
bool shader_t::set_uniform_color(int loc, colorRGBA const &val) { // same as vector3d
	if (loc >= 0) {glUniform4fv(loc, 1, &val.R); return 1;} else {return 0;}
}
bool shader_t::set_uniform_color(int loc, colorRGB  const &val) { // same as vector3d
	if (loc >= 0) {glUniform3fv(loc, 1, &val.R); return 1;} else {return 0;}
}
bool shader_t::set_uniform_matrix_3x3(int loc, float const *const m, bool transpose, unsigned num) {
	assert(num > 0);
	if (loc >= 0) {glUniformMatrix3fv(loc, num, transpose, m); return 1;} else {return 0;}
}
bool shader_t::set_uniform_matrix_4x4(int loc, float const *const m, bool transpose, unsigned num) {
	assert(num > 0);
	if (loc >= 0) {glUniformMatrix4fv(loc, num, transpose, m); return 1;} else {return 0;}
}

bool shader_t::add_uniform_float_array(char const *const name, float const *const val, unsigned num) const {
	return set_uniform_float_array(get_uniform_loc(name), val, num);
}
bool shader_t::add_uniform_float(char const *const name, float val) const {
	return set_uniform_float(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_int(char const *const name, int val) const {
	return set_uniform_int(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_uint(char const *const name, unsigned val) const {
	return set_uniform_uint(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_handle(char const *const name, GLuint64 val) const {
	return set_uniform_handle(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_vector2d(char const *const name, vector2d const &val) const {
	return set_uniform_vector2d(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_vector3d(char const *const name, vector3d const &val) const {
	return set_uniform_vector3d(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_vector4d(char const *const name, vector4d const &val) const {
	return set_uniform_vector4d(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_color(char const *const name, colorRGBA const &val) const {
	return set_uniform_color(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_color(char const *const name, colorRGB  const &val) const {
	return set_uniform_color(get_uniform_loc(name), val);
}
bool shader_t::add_uniform_matrix_4x4(char const *const name, float const *m, bool transpose, unsigned num) const {
	return set_uniform_matrix_4x4(get_uniform_loc(name), m, transpose, num);
}


// *** subroutines ***


int shader_t::subroutine_val_t::get_ix_for_name(char const *const name) const {

	for (auto i = name_to_ix.begin(); i != name_to_ix.end(); ++i) {
		if (i->first == name) return i->second;
	}
	return -1; // not found
}

unsigned shader_t::get_subroutine_index(int shader_type, char const *const name) const {

	assert(program && name);
	assert(shader_type < (int)NUM_SHADER_TYPES);
	assert(glGetSubroutineIndex != nullptr);
	unsigned const ret(glGetSubroutineIndex(program, shader_type_table[shader_type], name));
	assert(ret != GL_INVALID_INDEX); // check that name exists
	return ret;
}

unsigned shader_t::get_subroutine_uniform_loc(int shader_type, char const *const name) const {

	assert(program && name);
	assert(shader_type < (int)NUM_SHADER_TYPES);
	assert(glGetSubroutineUniformLocation != nullptr);
	unsigned const ret(glGetSubroutineUniformLocation(program, shader_type_table[shader_type], name));
	assert(ret != GL_INVALID_INDEX); // check that name exists
	return ret;
}

void shader_t::set_all_subroutines(int shader_type, unsigned count, char const *const *const uniforms, char const *const *const bindings) {

	assert(shader_type < (int)NUM_SHADER_TYPES && count > 0 && uniforms != nullptr && bindings != nullptr);
	if (glGetSubroutineIndex == nullptr) return;
	subroutine_val_t &val(subroutines[shader_type]);
	val.resize(count);
	
	for (unsigned i = 0; i < count; ++i) {
		unsigned const loc(get_subroutine_uniform_loc(shader_type, uniforms[i]));
		assert(loc < count);
		val.ixs[loc] = get_subroutine_index(shader_type, bindings[i]);
		int const ix(val.get_ix_for_name(uniforms[i]));
		assert(ix == -1); // no duplicate names
		val.name_to_ix.push_back(make_pair(uniforms[i], loc)); // add a new name=>ix mapping
	}
	set_subroutines(shader_type, val.ixs);
}

void shader_t::set_subroutines(int shader_type, unsigned count, unsigned const *const indices) {

	assert(program);
	assert(shader_type < (int)NUM_SHADER_TYPES && count > 0 && indices != nullptr);
	if (glGetSubroutineIndex == nullptr) return;
	int const stype(shader_type_table[shader_type]);
	//assert(count == glGetProgramStageiv(program, stype, GL_ACTIVE_SUBROUTINE_UNIFORM_LOCATIONS, v));
	//assert(all_indices < glGetProgramStageiv(program, stype, GL_ACTIVE_SUBROUTINES, v));
	assert(glUniformSubroutinesuiv != nullptr);
	glUniformSubroutinesuiv(stype, count, indices);
}

void shader_t::set_subroutine(int shader_type, unsigned index) { // single subroutine version

	subroutine_val_t &val(subroutines[shader_type]);
	val.resize(1); // must be empty or already size 1
	val.ixs[0] = index; // what about val.name_to_ix?
	set_subroutines(shader_type, val.ixs);
}

void shader_t::reset_subroutine(int shader_type, char const *const uniform, char const *const binding) {

	assert(uniform != nullptr && binding != nullptr);
	if (glGetSubroutineIndex == nullptr) return;
	subroutine_val_t &val(subroutines[shader_type]);
	int const ix(val.get_ix_for_name(uniform));
	assert(ix >= 0); // name must be known
	assert((unsigned)ix < val.ixs.size());
	val.ixs[ix] = get_subroutine_index(shader_type, binding); // cache this too?
	set_subroutines(shader_type, val.ixs);
}

void shader_t::restore_subroutines() {
	for (auto i = subroutines.begin(); i != subroutines.end(); ++i) {set_subroutines(i->first, i->second.ixs);}
}


// *** attrib variables setup ***


int shader_t::get_attrib_loc(char const *const name, bool allow_fail) const {

	assert(program && name);
	int const loc(glGetAttribLocation(program, name));
	if (!(allow_fail || loc >= 0)) {cerr << "Error: Failed to get shader attribute location '" << name << "'." << endl;}
	assert(allow_fail || loc >= 0); // Note: if variable is unused, loc will be -1
	return loc;
}


void shader_t::register_attrib_name(char const *const name, unsigned bind_ix) {

	assert(bind_ix < 100); // sanity check
	int const loc(get_attrib_loc(name)); // okay if -1
	if (bind_ix >= attrib_locs.size()) {attrib_locs.resize(bind_ix+1);}
	attrib_locs[bind_ix] = loc;
}


int shader_t::attrib_loc_by_ix(unsigned ix, bool allow_fail) const {

	if (allow_fail && ix >= attrib_locs.size()) return -1; // not set
	assert(ix < attrib_locs.size());
	return attrib_locs[ix]; // is it legal for this to return -1?
}


bool shader_t::set_attrib_float_array(int loc, float const *const val, unsigned num) const {

	if (loc < 0) {return 0;}
	switch (num) {
		case 1: glVertexAttrib1fv(loc, val); break;
		case 2: glVertexAttrib2fv(loc, val); break;
		case 3: glVertexAttrib3fv(loc, val); break;
		case 4: glVertexAttrib4fv(loc, val); break;
		case 16:
			for (unsigned n = 0; n < 4; ++n) {glVertexAttrib4fv(loc+n, val+4*n);}
			break; // mat4
		default: assert(0);
	}
	return 1;
}

bool shader_t::set_attrib_int_array(int loc, int const *const val, unsigned num) const { // Note: unused

	if (loc < 0) {return 0;}
	switch (num) {
	case 1: glVertexAttribI1iv(loc, val); break;
	case 2: glVertexAttribI2iv(loc, val); break;
	case 3: glVertexAttribI3iv(loc, val); break;
	case 4: glVertexAttribI4iv(loc, val); break;
	case 16:
		for (unsigned n = 0; n < 4; ++n) {glVertexAttribI4iv(loc+n, val+4*n);}
		break; // mat4
	default: assert(0);
	}
	return 1;
}

bool shader_t::set_attrib_float(int loc, float val) const {
	if (loc >= 0) {glVertexAttrib1f(loc, val); return 1;} else {return 0;}
}
bool shader_t::set_attrib_int(int loc, int val) const {
	if (loc >= 0) {glVertexAttrib1s(loc, val); return 1;} else {return 0;}
}

bool shader_t::add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const {
	return set_attrib_float_array(attrib_loc_by_ix(ix), val, num);
}
bool shader_t::add_attrib_int_array(unsigned ix, int const *const val, unsigned num) const {
	return set_attrib_int_array(attrib_loc_by_ix(ix), val, num);
}

bool shader_t::add_attrib_float(unsigned ix, float val) const {
	return set_attrib_float(attrib_loc_by_ix(ix), val);
}
bool shader_t::add_attrib_int(unsigned ix, int val) const {
	return set_attrib_int(attrib_loc_by_ix(ix), val);
}


// *** other variables setup ***


void shader_t::setup_enabled_lights(unsigned num, unsigned shaders_enabled) {

	assert(num <= MAX_SHADER_LIGHTS);
	if (prog_name_prefix.empty()) {prog_name_prefix.reserve(num+2);}
	prog_name_prefix.push_back(',');
	prog_name_prefix.push_back('L');
	char name_0[] = "const bool enable_light0 = false;";
	char name_1[] = "const bool enable_light0 = true;";

	for (unsigned i = 0; i < num; ++i) { // 0=sun, 1=moon, ...
		bool const enabled(is_light_enabled(i));
		prog_name_prefix.push_back(enabled ? '1' : '0');
		char *name(enabled ? name_1 : name_0);
		name[23] = char('0'+i);

		for (unsigned s = 0; s < NUM_SHADER_TYPES; ++s) { // put into correct shader(s): V, F, G, TC, TE, C
			if (shaders_enabled & (1<<s)) {set_prefix(name, s);}
		}
	}
}

void shader_t::upload_light_source(unsigned light_id, unsigned field_filt) {

	assert(is_setup());
	assert(light_id < MAX_SHADER_LIGHTS); // only supporting 8 light sources
	if (field_filt == 0 || !is_light_enabled(light_id)) return;
	light_loc_t &lloc(light_locs[light_id]);

	if (!lloc.valid) {
		static char ls_strs[MAX_SHADER_LIGHTS][5][40] = {};
		static bool is_setup(0);

		if (!is_setup) {
			for (unsigned i = 0; i < MAX_SHADER_LIGHTS; ++i) {
				sprintf(ls_strs[i][0], "fg_LightSource[%u].position", i);
				sprintf(ls_strs[i][1], "fg_LightSource[%u].ambient",  i);
				sprintf(ls_strs[i][2], "fg_LightSource[%u].diffuse",  i);
				sprintf(ls_strs[i][3], "fg_LightSource[%u].specular", i);
				sprintf(ls_strs[i][4], "fg_LightSource[%u].atten",    i);
			}
			is_setup = 1;
		}
		for (unsigned i = 0; i < 5; ++i) {lloc.v[i] = get_uniform_loc(ls_strs[light_id][i]);}
		lloc.valid = 1;
	}
	gl_light_params_t const &lp(gl_light_params[light_id]);
	gl_light_params_t &plp(prev_lps[light_id]);

	if ((field_filt & 0x01) && lloc.v[0] >= 0 && (lp.eye_space_pos != plp.eye_space_pos || lp.pos_w != plp.pos_w)) {
		vector4d const pos(lp.eye_space_pos, lp.pos_w);
		glUniform4fv(lloc.v[0], 1, &pos.x); // in eye space, using MVM at the time set_gl_light_pos() was called
		plp.eye_space_pos = lp.eye_space_pos; plp.pos_w = lp.pos_w;
	}
	if ((field_filt & 0x02) && lloc.v[1] >= 0 && lp.ambient  != plp.ambient ) {glUniform4fv(lloc.v[1], 1, &lp.ambient.R ); plp.ambient  = lp.ambient; }
	if ((field_filt & 0x04) && lloc.v[2] >= 0 && lp.diffuse  != plp.diffuse ) {glUniform4fv(lloc.v[2], 1, &lp.diffuse.R ); plp.diffuse  = lp.diffuse; }
	if ((field_filt & 0x08) && lloc.v[3] >= 0 && lp.specular != plp.specular) {glUniform4fv(lloc.v[3], 1, &lp.specular.R); plp.specular = lp.specular;}
	if ((field_filt & 0x10) && lloc.v[4] >= 0 && lp.atten    != plp.atten   ) {glUniform3fv(lloc.v[4], 1, &lp.atten.x   ); plp.atten    = lp.atten;   }
}

void shader_t::upload_light_sources_range(unsigned start, unsigned end) {

	for (unsigned i = start; i < end; ++i) {
		if (is_light_enabled(i)) {upload_light_source(i);}
	}
}

void shader_t::setup_scene_bounds() const {
	setup_scene_bounds_from_bcube(get_scene_bounds_bcube());
	add_uniform_vector3d("camera_pos", get_camera_pos());
}
void shader_t::setup_scene_bounds_from_bcube(cube_t const &bcube) const {
	add_uniform_vector3d("scene_llc",   bcube.get_llc());
	add_uniform_vector3d("scene_scale", bcube.get_size());
}

void shader_t::setup_fog_scale() const {
	add_uniform_float("fog_scale", (fog_enabled ? 1.0 : 0.0));
	add_uniform_float("fog_end",   cur_fog_end);
	add_uniform_color("fog_color", cur_fog_color);
}

void shader_t::check_for_fog_disabled() {
	if (!fog_enabled) {for (unsigned d = 0; d < 2; ++d) {set_prefix("#define NO_FOG", d);}} // VS/FS
}


void shader_t::set_prefix(char const *const prefix, unsigned shader_type) {

	assert(shader_type < NUM_SHADER_TYPES);
	if (prog_name_prefix.empty()) {prog_name_prefix.reserve(strlen(prefix) + 4);} // optimization
	prog_name_prefix.push_back(',');
	prog_name_prefix.push_back('s');
	prog_name_prefix.push_back('0'+shader_type);
	prog_name_prefix += prefix;
	prepend_string[shader_type] += prefix;
	prepend_string[shader_type].push_back('\n');
}

void shader_t::set_prefixes(char const *const prefix, unsigned shaders_enabled) {

	for (unsigned s = 0; s < NUM_SHADER_TYPES; ++s) { // put into correct shader(s): V, F, G, TC, TE, C
		if (shaders_enabled & (1<<s)) {set_prefix(prefix, s);}
	}
}

void shader_t::set_int_prefix(char const *const name, int val, unsigned shader_type) {
	
	if (val >= 0 && val <= 9) { // faster single character optimization
		set_prefix_str((string("const int ") + name + '=' + char('0' + val) + ';'), shader_type);
	}
	else {
		ostringstream oss;
		oss << val;
		set_prefix_str((string("const int ") + name + '=' + oss.str() + ';'), shader_type);
	}
}


void shader_t::set_color_e(colorRGBA const &color) {
	ensure_uniform_loc(emission_loc, "emission");
	set_uniform_color(emission_loc, color);
}

void shader_t::set_black_diffuse_emissive_color(colorRGBA const &color) {
	set_color_e(color);
	set_cur_color((emission_loc >= 0) ? colorRGBA(BLACK, color.alpha) : color); // if emissive color was set, then use black; otherwise, assume lighting is disabled and use the provided color
}

void shader_t::set_specular_color(colorRGB const &specular, float shininess) {

	assert(is_setup());
	shininess = max(0.0f, min(128.0f, shininess));
	colorRGBA spec_shine(specular, shininess); // pack specular.rgb and shininess together
	if (is_cloudy && world_mode != WMODE_UNIVERSE) {spec_shine *= 0.5;}

	if (spec_shine != last_spec) {
		ensure_uniform_loc(specular_color_loc, "specular_color");
		set_uniform_color(specular_color_loc, spec_shine);
		last_spec = spec_shine;
	}
}

void shader_t::set_material(base_mat_t const &mat) {
	set_specular_color(mat.spec_color, mat.shine);
	set_cur_color(mat.color);
}


// *** shader and program setup ***


struct program_t {
	unsigned p, sixs[NUM_SHADER_TYPES] = {0};
	bool valid;

	program_t() : p(0), valid(0) {}
	program_t(unsigned p_, unsigned sixs_[NUM_SHADER_TYPES]) : p(p_), valid(1) {}
};


class string_prog_map : public map<string, program_t> {
public:
	void clear() {
		free_data();
		map<string, program_t>::clear();
	}
	void free_data() {
		for (const_iterator s = begin(); s != end(); ++s) {glDeleteProgram(s->second.p);}
	}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~string_prog_map() {free_data();}
};


struct ix_valid_t {
	unsigned ix;
	bool valid;
	ix_valid_t() : ix(0), valid(0) {}
	ix_valid_t(unsigned const ix_) : ix(ix_), valid(1) {}
};


class string_shad_map : public map<string, ix_valid_t> {
public:
	void clear() {
		free_data();
		map<string, ix_valid_t>::clear();
	}
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			//assert(glIsShader(i->second.ix));
			glDeleteShader(i->second.ix);
		}
	}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~string_shad_map() {free_data();}
};


class shader_manager_t {

	string_prog_map loaded_programs;
	string_shad_map loaded_shaders[NUM_SHADER_TYPES]; // vertex=0, fragment=1, geometry=2, tess_control=3, tess_eval=4, compute=5
	map<string, string> loaded_files;
	string shader_id_str; // here to avoid constant reallocation

public:
	void clear() {
		loaded_programs.clear();
		for (unsigned d = 0; d < NUM_SHADER_TYPES; ++d) {loaded_shaders[d].clear();}
	}

	void clear_and_reload() {
		clear();
		loaded_files.clear();
	}

	bool check_strip_wrapper_chars(string &str, char cs, char ce) {
		if (str.front() == cs && str.back() == ce) {str = str.substr(1, str.size()-2); return 1;} // strip off the characters
		return 0;
	}
	bool load_shader_file(string const &fname, string &data, set<string> &all_fns) {
		if (fname.empty()) return 0;
		auto i(loaded_files.find(fname));
	
		if (i != loaded_files.end()) {
			data += i->second;
			return 1;
		}
		ifstream in(fname.c_str());
		if (!in.good()) return 0;
		string line, file_contents, str;
		
		while (std::getline(in, line)) {
			if (line.size() > 8) { // look for include directive
				istringstream iss(line);
				if ((iss >> str) && str == "#include") {
					// Note: it may help with MSVS syntax highlighting to add "#extension GL_ARB_shading_language_include : enable" to the shader,
					// though it likely doesn't understand the shaders directory system or ".part" rather than ".vert", ".frag" extensions
					// alternatively, there appears to be an ARB extension for shader includes:
					//glNamedStringARB(GL_SHADER_INCLUDE_ARB, -1, "header.glsl", -1, header_str);
					//glCompileShaderIncludeARB(shader, inc_dirs, _countof(inc_dirs), NULL);
					if (!(iss >> str)) {cerr << "Error: empty shader include" << endl; return 0;}
					// does this need to handle quoted string with spaces?
					if      (check_strip_wrapper_chars(str, '<',  '>' )) {str = "shaders/" + str;} // strip off the angle brackets (include from shader directory)
					else if (check_strip_wrapper_chars(str, '\"', '\"')) {} // strip off the quotes (include from local directory)
					if (str == fname) {cerr << "Error: recursive include of shader file '" << fname << "'" << endl; return 0;}
					// load the contents of this shader directly into the file contents of the including shader (inline it)
					if (!load_shader_file(str, file_contents, all_fns)) {cerr << "Error: Failed to load included shader file '" << str << "'" << endl; return 0;}
					file_contents += '\n';
					all_fns.insert(str); // record so that this file can be reloaded on error
					continue;
				}
			}
			file_contents += line + '\n';
		}
		loaded_files[fname] = file_contents;
		if (PRINT_SHADER) cout << "shader data:" << endl << file_contents << endl;
		data += file_contents;
		return 1;
	}

	bool clear_shader_file(string const &fname) {
		assert(!fname.empty());
		auto i(loaded_files.find(fname));
		if (i == loaded_files.end()) return 0; // not found - error?
		loaded_files.erase(i);
		return 1;
	}

	static void filename_split(string const &fname, vector<string> &fns, char sep) {
		stringstream ss(fname);
		string fn;
		while (getline(ss, fn, sep)) {fns.push_back(fn);}
	}

	static void get_shader_filenames(string const &name, unsigned type, vector<string> &fns) {
		assert(type < NUM_SHADER_TYPES);
		if (!shader_prefix_shared_file.empty()) {fns.push_back(shader_prefix_shared_file);} // add first, if nonempty
		if (!shader_prefix_files[type].empty()) {fns.push_back(shader_prefix_files[type]);} // add next, if nonempty
		filename_split(name, fns, '+');

		for (vector<string>::iterator i = fns.begin(); i != fns.end(); ++i) {
			assert(!i->empty());
			string fname(shaders_dir + "/" + *i);
		
			if ((*i)[i->size()-1] == '*') { // wildcard shader file: works with all shader types
				fname.erase(fname.size()-1);
			}
			else { // add shader type extension
				fname += "." + shader_name_table[type];
			}
			*i = fname;
		}
	}

	ix_valid_t &get_shader_by_name(string const &name, unsigned type) {return loaded_shaders[type][name];}
	program_t &get_program_by_name(string const &name) {return loaded_programs[name];}

	program_t &get_program_by_shader_names(string const shader_names[NUM_SHADER_TYPES], string const &prefix) {
		shader_id_str.clear();
		shader_id_str += prefix;
		for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {shader_id_str += shader_names[i] + ',';} // unique program identifier
		return get_program_by_name(shader_id_str);
	}
	void print_stats() const {
		cout << "Loaded shader files: " << loaded_files.size() << endl;
		cout << "Shader counts:";
		for (unsigned n = 0; n < NUM_SHADER_TYPES; ++n) {cout << " " << loaded_shaders[n].size();}
		cout << endl;
		cout << "Shader programs: " << loaded_programs.size() << endl;
	}
};

shader_manager_t shader_manager;


bool setup_shaders() {

	if (0) {
		int max_pv(0), max_level(0), def_o(0), def_i(0);
		glGetIntegerv(GL_MAX_PATCH_VERTICES, &max_pv);
		glGetIntegerv(GL_MAX_TESS_GEN_LEVEL, &max_level);
		glGetIntegerv(GL_PATCH_DEFAULT_OUTER_LEVEL, &def_o);
		glGetIntegerv(GL_PATCH_DEFAULT_INNER_LEVEL, &def_i);
		cout << "max_pv: " << max_pv << ", max_level: " << max_level << ", def levels: " << def_o << " " << def_i << endl; // 32 64 1 1
		//glPatchParameteri(GL_PATCH_VERTICES, v);
		//glPatchParameteri(GL_PATCH_DEFAULT_OUTER_LEVEL, v);
		//glPatchParameteri(GL_PATCH_DEFAULT_INNER_LEVEL, v);
	}
	cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;
	cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
	cout << "GLSL Shader Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader && GL_EXT_geometry_shader4) return 1;
	cerr << "Error setting up vertex and fragment GLSL shaders." << endl;
	return 0;
}


void clear_shaders() {
	clear_cached_shaders();
	shader_manager.clear();
}

void reload_all_shaders() { // clears and reloads *everything*
	// Note: do we want/need some function called every frame that check if shader files have been modified and calls this?
	cout << "Reloading all shaders" << endl;
	clear_cached_shaders();
	shader_manager.clear_and_reload();
	clear_tiled_terrain_shaders();
}

void print_shader_stats() {shader_manager.print_stats();}


bool yes_no_query(string const &query_str) {

	while (1) {
		cout << query_str << " [Y|N]" << endl;
		string str;
		cin >> str;
		if (str == "Y" || str == "y" || str == "1") {return 1;}
		if (str == "N" || str == "n" || str == "0") {return 0;}
		cout << "Answer not understood" << endl;
	}
	return 0; // never gets here
}


void write_shader_file(string const &name, string const &data, int type) {

	string const out_dir("shaders/generated");
	string filename;

	for (auto i = name.begin(); i != name.end(); ++i) {
		if (*i != '*') {filename.push_back(*i);} // skip special character
		if (filename.size() > 160) break; // limit filename length
	}
	string const out_fn(out_dir + "/" + filename + "." + shader_name_table[type]);
	cout << "writing shader file " << out_fn << endl;
	ofstream out(out_fn);
	assert(out.good());
	out << data << endl;
}


unsigned shader_t::get_shader(string const &name, unsigned type) const {
	
	//RESET_TIME;
	if (name.empty()) return 0; // none selected
	assert(type < NUM_SHADER_TYPES);
	string const lookup_name(name + prepend_string[type]);
	ix_valid_t &ixv(shader_manager.get_shader_by_name(lookup_name, type));
	if (ixv.valid) {return ixv.ix;} // already loaded
	
	// create a new shader
	string const version_info("#version 430\n"); // use version 430 for all shaders
	vector<string> fns;
	shader_manager.get_shader_filenames(name, type, fns);
	set<string> all_fns(fns.begin(), fns.end()); // fns + include files
	bool failed(0);
	unsigned shader(0);

	while (1) { // retry loop
		if (failed) {
			if (fullscreen || !yes_no_query("Retry?")) { // don't query the user when maximized
				cerr << "Exiting." << endl;
				exit(1);
			}
			for (set<string>::const_iterator i = all_fns.begin(); i != all_fns.end(); ++i) {
				if (shader_manager.clear_shader_file(*i)) {cout << "Reloading shader component " << *i << endl;}
			}
			failed = 0;
		}
		string data(version_info + prepend_string[type]);
		//data += "#line 1\n"; // add this to start from line 0 here to exclude header lines

		for (vector<string>::const_iterator i = fns.begin(); i != fns.end(); ++i) {
			if (!shader_manager.load_shader_file(*i, data, all_fns)) {
				cerr << "Error loading shader file " << *i << "." << endl;
				failed = 1; break;
			}
		}
		if (failed) continue;
		if (DEBUG_SHADER) {cout << "final shader data for <" << name << ">:" << endl << data << endl;}
		if (GEN_FINAL_SHADER_FILES) {write_shader_file(name, data, type);}
		assert(shader == 0);
		shader = glCreateShader(shader_type_table[type]);
	
		if (!shader) {
			cerr << "Error: Failed to create " << shader_name_table[type] << " shader " << name << "." << endl;
			failed = 1; continue;
		}
		//RESET_TIME;
		const char *src(data.c_str());
		glShaderSource(shader, 1, &src, 0);
		glCompileShader(shader);
		int status(0);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
		//PRINT_TIME("Compile Shader");

		if (status != GL_TRUE) {
			write_shader_file(name, data, type); // write it out so that we can reference the correct line numbers
			cerr << "Compilation of " << shader_name_table[type] << " shader " << name << " failed with status " << status << endl;
			print_shader_info_log(shader);
			cerr << endl;
			glDeleteShader(shader); shader = 0; // recreate rather than reusing to ensure there is no bad state (probably unnecessary)
			failed = 1; continue;
		}
		break;
	} // while(1)
	if (PRINT_LOG) {print_shader_info_log(shader);}
	ixv = ix_valid_t(shader); // cache the shader
	//PRINT_TIME("Create Shader");
	return shader;
}


// See http://www.lighthouse3d.com/tutorials/glsl-core-tutorial/glsl-core-tutorial-create-a-program/
bool shader_t::begin_shader(bool do_enable) {

	//RESET_TIME;
	// get the program
	program_t &prog(shader_manager.get_program_by_shader_names(shader_names, prog_name_prefix));

	if (prog.valid) { // program already exists
		program = prog.p;
	}
	else { // create a new program
		unsigned shader_ixs[NUM_SHADER_TYPES] = {0};
		
		while (1) { // retry loop
			check_gl_error(299);
			program = glCreateProgram();

			for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
				// Note: we *can* attach multiple shaders of the same type to a single program as long as only one has a main()
				shader_ixs[i] = get_shader(shader_names[i], i);
				if (shader_ixs[i]) {glAttachShader(program, shader_ixs[i]);}
			}
			// vertex and fragment shaders are required, geometry shader is optional; or only a compute shader
			assert((shader_ixs[0] && shader_ixs[1]) || shader_ixs[5]);
			if (shader_ixs[2]) {assert(GL_EXT_geometry_shader4);} // geometry shader
			check_gl_error(300);
			glLinkProgram(program);
			int status(0);
			glGetProgramiv(program, GL_LINK_STATUS, &status);
			check_gl_error(298);
			if (status == GL_TRUE) break; // success
			cerr << "Linking of program ";
			for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {cout << shader_names[i] + " ";}
			cerr << "failed with status " << status << endl;
			print_program_info_log();
			cerr << endl;
				
			if (!yes_no_query("Retry?")) {
				cerr << "Exiting." << endl;
				exit(1);
			}
			// do something simple and reload all shaders, rather than trying to figure out exactly what happened and what needs to be reloaded
			// somewhat primitive, but this should rarely happen, and gets the job done
			glDeleteProgram(program); // delete and recreate
			shader_manager.clear_and_reload();
		} // end retry loop
		if (PRINT_LOG) {print_program_info_log();}
#if 0 // debugging
		GLenum binary_format(0);
		vector<unsigned char> binary_data;
		get_program_binary(binary_data, binary_format);
		set_program_binary(binary_data, binary_format);
#endif
		for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
			if (shader_ixs[i]) {glDetachShader(program, shader_ixs[i]);}
		}
		prog = program_t(program, shader_ixs); // cache the program
		//PRINT_TIME("Create Program");
	}
	cache_vnct_locs();
	cache_matrix_locs();
	emission_loc = specular_color_loc = -1;
	if (do_enable) {enable();}
#if 0 // debugging
	glValidateProgram(program);
	int status(0);
	glGetProgramiv(program,  GL_VALIDATE_STATUS, &status);
	if (status == GL_FALSE) {cerr << "Error: Program validation failed" << endl;}
#endif
	return 1;
}

void shader_t::get_program_binary(vector<unsigned char> &binary_data, GLenum &binary_format) const {
	int prog_length(0), length(0);
	glGetProgramiv(program, GL_PROGRAM_BINARY_LENGTH, &prog_length);
	//cout << TXT(prog_length) << endl;
	assert(prog_length > 0);
	binary_data.resize(prog_length);
	glGetProgramBinary(program, prog_length, &length, &binary_format, &binary_data.front());
	assert(length == prog_length);
}
void shader_t::set_program_binary(vector<unsigned char> const &binary_data, GLenum const binary_format) {
	glProgramBinary(program, binary_format, &binary_data.front(), binary_data.size());
}


void shader_t::print_shader_info_log(unsigned shader) {

	int len(0), len2(0);
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);

	if (len > 0) {
		vector<char> info_log_msg(len);
		glGetShaderInfoLog(shader, len, &len2, &info_log_msg.front()); 
		assert(len2 <= len);
		cout << "Info log: " << string(info_log_msg.begin(), info_log_msg.end());
	}
}

void shader_t::print_program_info_log() const {

	int len(0), len2(0);
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);

	if (len > 0) {
		vector<char> info_log_msg(len);
		glGetProgramInfoLog(program, len, &len2, &info_log_msg.front()); 
		assert(len2 <= len);
		cout << "Info log: " << string(info_log_msg.begin(), info_log_msg.end());
	}
}


void shader_t::end_shader() { // ok to call if not in a shader
	disable();
	clear();
}

void shader_t::clear() {

	program = 0;
	
	for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
		prepend_string[i].clear();
		shader_names[i].clear();
	}
	for (unsigned i = 0; i < MAX_SHADER_LIGHTS; ++i) {
		light_locs[i].invalidate();
		prev_lps[i] = gl_light_params_t(); // reset
	}
	prog_name_prefix.clear();
	attrib_locs.clear();
	subroutines.clear();
	clear_vntc_locs();
	clear_properties();
	last_spec  = ALPHA0;
	user_flags = 0; // clear user flags
}


void shader_t::make_current() {
	assert(program);
	glUseProgram(program);
	restore_subroutines();
	cur_shader = this;
}

void shader_t::enable() {
	make_current();
	upload_all_light_sources();
	upload_pjm();
	upload_mvm();
	mvm_changed = 0;
}
void shader_t::disable() {
	cur_shader = NULL; // must be done first to prevent check_mvm_update() from trying to upload a new MVM in enable_vnct_atribs()
	if (is_setup()) {enable_vnct_atribs(0, 0, 0, 0);} // disable all
	glUseProgram(0);
}


bool shader_is_active() {return (cur_shader != nullptr);}

// Note 1: We don't handle the case where the projection matrix is updated while a shader is active because this currently doesn't occur.
// Note 2: We assume we only need to update the MVM in at most one active shader, and once updated we can clear the changed flag
void check_mvm_update() {
	if (!mvm_changed) return; // nothing to update
	if (cur_shader) {cur_shader->upload_mvm();}
	mvm_changed = 0;
}

void shader_t::cache_matrix_locs() {
	pm_loc   = get_uniform_loc("fg_ProjectionMatrix"); // okay if returns -1
	mvm_loc  = get_uniform_loc("fg_ModelViewMatrix");
	mvmi_loc = get_uniform_loc("fg_ModelViewMatrixInverse");
	mvpm_loc = get_uniform_loc("fg_ModelViewProjectionMatrix");
	nm_loc   = get_uniform_loc("fg_NormalMatrix");
}

void shader_t::upload_pjm() { // projection matrix
	if (pm_loc >= 0) {set_uniform_matrix_4x4(pm_loc, fgGetPJM().get_ptr(), 0);} // transpose = 0
}

xform_matrix xform_matrix::inverse() const {return glm::affineInverse((glm::mat4)*this);}

void shader_t::upload_mvm() { // and everything that depends on the mvm

	if (mvm_loc < 0 && mvmi_loc < 0 && mvpm_loc < 0 && nm_loc < 0) return; // nothing to update
	xform_matrix const &mvm(fgGetMVM());
	if (mvm_loc >= 0) {set_uniform_matrix_4x4(mvm_loc, mvm.get_ptr(), 0);} // modelview matrix

	if (mvmi_loc >= 0) { // inverse modelview matrix
		set_uniform_matrix_4x4(mvmi_loc, mvm.inverse().get_ptr(), 0);
	}
	if (mvpm_loc >= 0) { // modelview-projection matrix
		xform_matrix const mvp(fgGetPJM() * mvm);
		set_uniform_matrix_4x4(mvpm_loc, mvp.get_ptr(), 0);
	}
	if (nm_loc >= 0) { // normal matrix
		glm::mat3 const nm(glm::inverseTranspose(glm::mat3(mvm)));
		set_uniform_matrix_3x3(nm_loc, glm::value_ptr(nm), 0);
	}
}


// built-in attribute setup/enable/binding
void shader_t::cache_vnct_locs() { // Note: program need not be enabled
	assert(is_setup());
	// Note: locations will generally be {0,1,2,3} as assigned in the vertex shader, but if unused they can be -1
	const char *loc_strs[4] = {"fg_Vertex", "fg_Normal", "fg_Color", "fg_TexCoord"};
	for (unsigned i = 0; i < 4; ++i) {vnct_locs[i] = get_attrib_loc(loc_strs[i], 1);} // okay if fails
}

void shader_t::enable_vnct_atribs(bool va, bool tca, bool na, bool ca) const { // Note: program must be enabled

	assert(is_setup());
	bool const enables[4] = {va, na, ca, tca};

	for (unsigned i = 0; i < 4; ++i) {
		int const loc(vnct_locs[i]);
		if (loc < 0) continue; // unused in this shader
		if (enables[i]) {glEnableVertexAttribArray(loc);} else {glDisableVertexAttribArray(loc);}
	}
	check_mvm_update();
}

void shader_t::set_vertex_ptr(unsigned stride, void const *const ptr) const {
	assert(vnct_locs[0] >= 0); // vertex must always be available
	glVertexAttribPointer(vnct_locs[0], 3, GL_FLOAT, GL_FALSE, stride, ptr);
}

void shader_t::set_normal_ptr(unsigned stride, void const *const ptr, bool compressed) const {
	if (vnct_locs[1] >= 0) {glVertexAttribPointer(vnct_locs[1], 3, (compressed ? GL_BYTE          : GL_FLOAT), compressed, stride, ptr);}
}

void shader_t::set_color4_ptr(unsigned stride, void const *const ptr, bool compressed) const {
	if (vnct_locs[2] >= 0) {glVertexAttribPointer(vnct_locs[2], 4, (compressed ? GL_UNSIGNED_BYTE : GL_FLOAT), compressed, stride, ptr);}
}

void shader_t::set_tcoord_ptr(unsigned stride, void const *const ptr, bool compressed) const {
	if (vnct_locs[3] >= 0) {glVertexAttribPointer(vnct_locs[3], 2, (compressed ? GL_SHORT         : GL_FLOAT), compressed, stride, ptr);}
}

void shader_t::set_cur_color(colorRGBA const &color) const {
	if (vnct_locs[2] >= 0) {glVertexAttrib4fv(vnct_locs[2], &color.R);}
}
void shader_t::set_cur_normal(vector3d const &normal) const {
	if (vnct_locs[1] >= 0) {glVertexAttrib3fv(vnct_locs[1], &normal.x);}
}


// some simple shared shaders
void shader_t::begin_color_only_shader() {
	set_vert_shader("vert_xform_only");
	set_frag_shader("color_only");
	begin_shader();
}
void shader_t::begin_color_only_shader(colorRGBA const &color) {
	begin_color_only_shader();
	set_cur_color(color);
}

bool is_csm_active();
void shader_csm_render_setup(shader_t &s);

void shader_t::begin_shadow_map_shader(bool use_alpha_mask, bool enable_xlate_scale) {
	bool const use_csm(is_csm_active());

	if (use_alpha_mask) {
		if (use_csm) {
			set_vert_shader("shadow_map_csm_tc");
			set_geom_shader("csm_layers_tc");
		}
		else {
			set_vert_shader("pos_only");
		}
		set_frag_shader("alpha_mask_shadow");
		begin_shader();
		add_uniform_float("min_alpha", MIN_SHADOW_ALPHA);
		add_uniform_int("tex0", 0);
	}
	else {
		if (use_csm) {
			set_vert_shader(enable_xlate_scale ? "shadow_map_csm_xlate_scale" : "shadow_map_csm");
			set_geom_shader("csm_layers");
		}
		else {
			set_vert_shader(enable_xlate_scale ? "vertex_xlate_scale" : "shadow_map");
		}
		set_frag_shader("empty_shader");
		begin_shader();
	}
	if (use_csm) {shader_csm_render_setup(*this);}
}

void shader_t::begin_simple_textured_shader(float min_alpha, bool include_2_lights, bool use_texgen, colorRGBA const *const color) {

	if (include_2_lights) {
		setup_enabled_lights(2, 1); // sun and moon VS lighting
		set_vert_shader(use_texgen ? "ads_lighting.part*+texture_gen.part+two_lights_texture_gen" : "ads_lighting.part*+two_lights_texture");
	}
	else {
		set_vert_shader(use_texgen ? "texture_gen.part+no_lighting_texture_gen" : "no_lighting_tex_coord");
	}
	bool const use_fog(include_2_lights && fog_enabled);
	set_frag_shader(use_fog ? "linear_fog.part+textured_with_fog" : "simple_texture");
	begin_shader();
	if (color) {set_cur_color(*color);}
	add_uniform_float("min_alpha", min_alpha);
	add_uniform_int("tex0", 0);
	if (use_fog) {setup_fog_scale();}
}

void shader_t::begin_untextured_lit_glcolor_shader() {
	begin_simple_textured_shader(0.0, 1); // lighting (not actually textured)
	select_texture(WHITE_TEX); // untextured
}


// compute shader


bool compute_shader_base_t::setup_target_texture(unsigned &tid, bool is_R32F) const {
	if (tid > 0) return 0; // already setup
	setup_texture(tid, 0, 0, 0, 0, 0, 1); // nearest, clamp, no mipmaps

	if (is_R32F) {
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, xsize, ysize, 0, GL_RED, GL_FLOAT, nullptr);
	}
	else {
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, xsize, ysize, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	}
	return 1;
}

// factored out so that we can make it virtual and override it in the future
void compute_shader_t::draw_geom() const {
	float const z = 0.0;
	bind_vao(0); // bind default VBO in core_context mode - required since this is called outside the display() loop
	vert_wrap_t verts[4] = {point(0,0,z), point(0,1,z), point(1,1,z), point(1,0,z)};
	draw_verts(verts, 4, GL_TRIANGLE_FAN); // one quad from [0,0] to [1,1] that exactly covers the viewport
}

void compute_shader_t::begin() {
	set_vert_shader("vert_xform_vpos");
	set_frag_shader(frag_shader_str);
	begin_shader();
}
void compute_shader_t::end_shader() {
	shader_t::end_shader();
	free_fbo(fbo_id);
}

void compute_shader_t::unset_fbo(bool keep_fbo_for_reuse) { // call once after run() calls
	if (!keep_fbo_for_reuse) {free_fbo(fbo_id);}
	disable_fbo();
}

void compute_shader_t::setup_and_run(unsigned &tid, bool is_R32F, bool is_first, bool is_last) {
	
	setup_target_texture(tid, is_R32F);
	
	if (is_first) {
		glViewport(0, 0, xsize, ysize);
		fgMatrixMode(FG_PROJECTION);
		fgPushIdentityMatrix();
		fgOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0); // [0,0] to [1,1]
		fgMatrixMode(FG_MODELVIEW);
		fgPushIdentityMatrix();
	}
	assert(!is_running);
	run(tid);
	is_running = 1;

	if (is_last) {
		restore_prev_mvm_pjm_state();
		set_standard_viewport(); // restore state
	}
}

void compute_shader_t::run(unsigned &tid) { // call N times between pre_run() and post_run() calls
	// setup fbo and draw geometry
	enable_fbo(fbo_id, tid, 0);
	set_temp_clear_color(BLACK);
	ensure_filled_polygons();
	draw_geom();
	reset_fill_mode();
}

// tid may or may not be setup prior to this call
void compute_shader_t::gen_matrix_RGBA8(vector<float> &vals, unsigned &tid, bool is_first, bool is_last, bool keep_fbo_for_reuse) {
	
	setup_and_run(tid, 0, is_first, is_last);
	vals.resize(xsize*ysize);
	vector<unsigned char> data(4*vals.size());
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, xsize, ysize, GL_RGBA, GL_UNSIGNED_BYTE, &data.front()); // GL_BGRA?

	for (unsigned y = 0; y < ysize; ++y) {
		for (unsigned x = 0; x < xsize; ++x) {
			unsigned const ix(y*xsize + x);
			vals[ix] = data[ix<<2]; // red component [0.0, 1.0)
		}
	}
	is_running = 0;
	if (is_last) {unset_fbo(keep_fbo_for_reuse);}
}

// tid may or may not be setup prior to this call
void compute_shader_t::gen_matrix_R32F(vector<float> &vals, unsigned &tid, bool is_first, bool is_last, bool keep_fbo_for_reuse) {
	setup_and_run(tid, 1, is_first, is_last);
	prep_for_read_pixels(is_first);
	read_float_vals(vals, is_last, keep_fbo_for_reuse);
}

void compute_shader_t::read_float_vals(vector<float> &vals, bool is_last, bool keep_fbo_for_reuse) {

	bind_fbo(fbo_id);
	assert(is_running);
	vals.resize(xsize*ysize);
	read_pixels(vals, is_last); // Note: slower on old cards, faster on new ones
	//glReadBuffer(GL_COLOR_ATTACHMENT0); glReadPixels(0, 0, xsize, ysize, GL_RED, GL_FLOAT, &vals.front());
	is_running = 0;
	if (is_last) {unset_fbo(keep_fbo_for_reuse);}
}

void compute_shader_t::prep_for_read_pixels(bool is_first) {

	bind_fbo(fbo_id);

	if (is_first) {
		glGenBuffers(1, &pbo);
		bind_pbo(pbo);
		glBufferData(GL_PIXEL_PACK_BUFFER, get_pbo_size(), NULL, GL_STREAM_READ);
	}
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, xsize, ysize, GL_RED, GL_FLOAT, nullptr);
	bind_pbo(0);
	disable_fbo();
}

void compute_shader_t::read_pixels(vector<float> &vals, bool is_last) {

	assert(pbo); bind_pbo(pbo);
	void *ptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, get_pbo_size(), GL_MAP_READ_BIT); // Note: blocks until data is ready
	assert(ptr != nullptr);
	memcpy((void *)vals.data(), ptr, get_pbo_size());
	glUnmapBuffer(GL_PIXEL_PACK_BUFFER);

	if (is_last) {
		bind_pbo(0);
		glDeleteBuffers(1, &pbo);
		pbo = 0;
	}
}


void pad_to_block_sz(unsigned &sz, unsigned block_sz) {
	assert(block_sz > 0);
	if (sz % block_sz) {sz += block_sz - (sz % block_sz);}
}
unsigned div_by_block_sz(unsigned sz, unsigned block_sz) {
	//if (sz == 1) return 1;
	assert(block_sz > 0);
	assert((sz % block_sz) == 0);
	return sz/block_sz;
}

compute_shader_comp_t::compute_shader_comp_t(string const &cstr, unsigned xsize_, unsigned ysize_, unsigned zsize_, unsigned bsx, unsigned bsy, unsigned bsz) :
	compute_shader_base_t(xsize_, ysize_), zsize(zsize_), zsize_req(zsize), block_sz_x(bsx), block_sz_y(bsy), block_sz_z(bsz), comp_shader_str(cstr)
{
	assert(zsize > 0 && block_sz_x > 0 && block_sz_y > 0 && block_sz_z > 0);
	pad_to_block_sz(xsize, block_sz_x);
	pad_to_block_sz(ysize, block_sz_y);
	pad_to_block_sz(zsize, block_sz_z);
}

void compute_shader_comp_t::begin() {
	ostringstream oss;
	oss << "layout (local_size_x = " << block_sz_x << ", local_size_y = " << block_sz_y << ", local_size_z = " << block_sz_z << ") in;";
	set_prefix(oss.str().c_str(), 5);
	set_comp_shader(comp_shader_str);
	begin_shader();
	add_uniform_int("dest_tex", 0);
}

void compute_shader_comp_t::setup_and_run(unsigned &tid, bool is_R32F, bool is_first, bool is_last) { // Note: is_first and is_last are unused

	assert(is_R32F); // only supported mode
	check_gl_error(500);

	if (is_3d()) {
		if (tid == 0) {
			setup_3d_texture(tid, GL_NEAREST, GL_CLAMP_TO_EDGE);
			glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, xsize, ysize, zsize, 0, GL_RED, GL_FLOAT, nullptr);
		}
		bind_3d_texture(tid);
	}
	else {
		setup_target_texture(tid, 1);
		bind_2d_texture(tid);
	}
	glBindImageTexture(0, tid, 0, (is_3d() ? GL_TRUE : GL_FALSE), 0, GL_WRITE_ONLY, GL_R32F); // use image unit 0, layered if 3D texture
	// num threads = div_by_block_sz(xsize, block_sz_x)*div_by_block_sz(ysize, block_sz_y)*div_by_block_sz(zsize, block_sz_z)
	// block size  = block_sz_x*block_sz_y*block_sz_z
	glDispatchCompute(div_by_block_sz(xsize, block_sz_x), div_by_block_sz(ysize, block_sz_y), div_by_block_sz(zsize, block_sz_z));
	assert(!is_running);
	is_running = 1;
}

void compute_shader_comp_t::read_float_vals(vector<float> &vals, bool is_last, bool keep_fbo_for_reuse) { // Note: is_last and keep_fbo_for_reuse are unused

	assert(is_running);
	is_running = 0;
	vals.resize(xsize*ysize*zsize);
	glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
	glGetTexImage((is_3d() ? GL_TEXTURE_3D : GL_TEXTURE_2D), 0, GL_RED, GL_FLOAT, &vals.front());

	if (xsize != xsize_req || ysize != ysize_req || zsize != zsize_req) { // need to copy a smaller sub-range from a larger texture
		vector<float> temp(xsize_req*ysize_req*zsize_req);

		for (unsigned y = 0; y < ysize_req; ++y) {
			for (unsigned x = 0; x < xsize_req; ++x) {
				unsigned const temp_off((x + y*xsize_req)*zsize_req), vals_off((x + y*xsize)*zsize);
				for (unsigned z = 0; z < zsize_req; ++z) {temp[z + temp_off] = vals[z + vals_off];}
			}
		}
		vals.swap(temp);
	}
	check_gl_error(501);
}

void compute_shader_comp_t::gen_matrix_R32F(vector<float> &vals, unsigned &tid, bool is_first, bool is_last) {
	setup_and_run(tid, 1, is_first, is_last);
	read_float_vals(vals);
}


// **************** INSTANCING + TRANSFORMS ****************


void upload_mvm_to_shader(shader_t &s, char const *const var_name) {
	s.add_uniform_matrix_4x4(var_name, fgGetMVM().get_ptr(), 0);
}


// Note: assumes cur_vbo is currently bound by the caller, and will leave it bound after the call; indices can be NULL
void instance_render_t::draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo, int index_type, void const *const indices, unsigned first, unsigned cur_vao) {

	if (inst_xforms.empty()) return;
	assert(loc >= 0); // Note: could handle this case
	void const *vbo_ptr(get_dynamic_vbo_ptr(inst_xforms.front().get_ptr(), inst_xforms.size()*sizeof(xform_matrix)));
	shader_float_matrix_uploader<4,4>::enable(loc, 1, (float const *)vbo_ptr); // use hardware instancing
	if (cur_vbo) {bind_vbo(cur_vbo);}
	if (cur_vao) {bind_vao(cur_vao);}
	
	if (index_type != GL_NONE) { // indexed
		glDrawElementsInstanced(prim_type, count, index_type, indices, inst_xforms.size());
	}
	else {
		assert(indices == nullptr);
		glDrawArraysInstanced(prim_type, first, count, inst_xforms.size());
	}
	++num_frame_draw_calls;
	shader_float_matrix_uploader<4,4>::disable(loc);
	inst_xforms.clear();
}

/*static*/ template<unsigned M, unsigned N> void shader_float_matrix_uploader<M, N>::enable(int start_loc, int divisor, float const *const data) {
	assert(start_loc >= 0 && divisor >= 0);
  
	for (unsigned i = 0; i < N; ++i) {
		int const loc(start_loc + i);
		glEnableVertexAttribArray(loc);
		glVertexAttribPointer(loc, M, GL_FLOAT, GL_FALSE, M*N*sizeof(float), (const void *)(data + M*i));
		glVertexAttribDivisor(loc, divisor);
	}
}
/*static*/ template<unsigned M, unsigned N> void shader_float_matrix_uploader<M, N>::disable(int start_loc) {
	for (unsigned i = 0; i < N; ++i) {disable_instancing_for_shader_loc(start_loc + i);}
}

// explicit instantiations
template struct shader_float_matrix_uploader<4, 4>;
template struct shader_float_matrix_uploader<3, 1>;


void set_point_sprite_mode(bool enabled) { // Note: to be removed when using a core profile

	if (enabled) {
		if (!init_core_context) {glEnable(GL_POINT_SPRITE);}
		glEnable(GL_PROGRAM_POINT_SIZE);
	}
	else {
		if (!init_core_context) {glDisable(GL_POINT_SPRITE);}
		glDisable(GL_PROGRAM_POINT_SIZE);
	}
}

