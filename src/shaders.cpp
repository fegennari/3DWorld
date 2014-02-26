// 3D World - Vertex and Fragment GLSL Shader Framework
// by Frank Gennari
// 11/1/10
#include "shaders.h"
#include "mesh.h" // for scene bounds
#include "gl_ext_arb.h"
#include <fstream>

using namespace std;

bool const PRINT_SHADER = 0; // print shaders loaded from files
bool const DEBUG_SHADER = 0; // print final generated shaders
bool const PRINT_LOG    = 0;

string const shaders_dir = "shaders";

extern bool fog_enabled;
extern unsigned enabled_lights;


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


bool shader_t::set_uniform_float_array(int loc, float const *const val, unsigned num) {
	if (loc >= 0) {glUniform1fv(loc, num, val); return 1;} else {return 0;}
}

bool shader_t::set_uniform_float(int loc, float val) {
	if (loc >= 0) {glUniform1f(loc, val); return 1;} else {return 0;}
}

bool shader_t::set_uniform_int(int loc, int val) {
	if (loc >= 0) {glUniform1i(loc, val); return 1;} else {return 0;}
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

bool shader_t::set_uniform_matrid_4x4(int loc, float *m, bool transpose) {
	if (loc >= 0) {glUniformMatrix4fv(loc, 1, transpose, m); return 1;} else {return 0;}
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

bool shader_t::add_uniform_vector2d(char const *const name, vector2d const &val) const {
	return set_uniform_vector2d(get_uniform_loc(name), val);
}

bool shader_t::add_uniform_vector3d(char const *const name, vector3d const &val) const {
	return set_uniform_vector3d(get_uniform_loc(name), val);
}

bool shader_t::add_uniform_color(char const *const name, colorRGBA const &val) const {
	return set_uniform_color(get_uniform_loc(name), val);
}

bool shader_t::add_uniform_color(char const *const name, colorRGB  const &val) const {
	return set_uniform_color(get_uniform_loc(name), val);
}

bool shader_t::add_uniform_matrid_4x4(char const *const name, float *m, bool transpose) const {
	return set_uniform_matrid_4x4(get_uniform_loc(name), m, transpose);
}


// unused, unfinished
bool shader_t::set_uniform_buffer_data(char const *name, float const *data, unsigned size, unsigned &buffer_id) const {

	assert(program && name);
	assert(data && size);

	// There's only one uniform block.
	int const uniformBlockIndex(glGetUniformBlockIndex(program, "uniform_block"));
	if (uniformBlockIndex < 0) return 0;
	//cout << "ix: " << uniformBlockIndex << endl;

	// Associate the uniform block to binding point 0
	glUniformBlockBinding(program, uniformBlockIndex, 0);

	// Get the uniform block's size
	int uniformBlockSize;
	glGetActiveUniformBlockiv(program, uniformBlockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &uniformBlockSize);
	//cout << "block_size: " << uniformBlockSize << endl;

	// The data uniform will need to be updated, so we'll query its offset/size.
	unsigned index;
	int offset, buf_size;

	// First, get the index for the uniform
	glGetUniformIndices(program, 1, &name, &index);
	assert((int)index >= 0);

	// Use the index to query offset and size
	glGetActiveUniformsiv(program, 1, &index, GL_UNIFORM_OFFSET, &offset);
	glGetActiveUniformsiv(program, 1, &index, GL_UNIFORM_SIZE,   &buf_size);
	//cout << "index: " << index << ", offset: " << offset << ", size: " << buf_size << endl;
	assert(size <= (unsigned)buf_size);
	size = min(size, (unsigned)buf_size);

	// Create UBO
	if (buffer_id == 0 || !glIsBuffer(buffer_id)) glGenBuffers(1, &buffer_id);
	assert(buffer_id > 0);
	glBindBuffer(GL_UNIFORM_BUFFER, buffer_id);
	
	// We can use BufferData to upload our data to the shader, since we know it's in the std140 layout
	glBufferData(GL_UNIFORM_BUFFER, uniformBlockSize, NULL, GL_DYNAMIC_DRAW);

	// Bind constants to UBO binding point 0
	glBindBufferBase(GL_UNIFORM_BUFFER, 0, buffer_id);
	glBufferSubData(GL_UNIFORM_BUFFER, offset, size, data);
	glBindBuffer(GL_UNIFORM_BUFFER, 0); // unbind
	return 1;
}


// *** attrib variables setup ***


int shader_t::get_attrib_loc(char const *const name, bool allow_fail) const {

	assert(program && name);
	int const loc(glGetAttribLocation(program, name));
	assert(allow_fail || loc >= 0); // Note: if variable is unused, loc will be -1
	return loc;
}


void shader_t::register_attrib_name(char const *const name, unsigned bind_ix) {

	assert(bind_ix < 100); // sanity check
	int const loc(get_attrib_loc(name));
	if (bind_ix >= attrib_locs.size()) attrib_locs.resize(bind_ix+1);
	attrib_locs[bind_ix] = loc;
}


int shader_t::attrib_loc_by_ix(unsigned ix) const {

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
			for (unsigned n = 0; n < 4; ++n) {
				glVertexAttrib4fv(loc+n, val+4*n);
			}
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

bool shader_t::add_attrib_float(unsigned ix, float val) const {
	return set_attrib_float(attrib_loc_by_ix(ix), val);
}

bool shader_t::add_attrib_int(unsigned ix, int val) const {
	return set_attrib_int(attrib_loc_by_ix(ix), val);
}


// *** other variables setup ***


void shader_t::setup_enabled_lights(unsigned num, unsigned shaders_enabled) {

	assert(num <= 8);
	prog_name_suffix.push_back(',');
	prog_name_suffix.push_back('L');
	char name[14] = "enable_light0";

	for (unsigned i = 0; i < num; ++i) { // 0=sun, 1=moon, ...
		bool const enabled(is_light_enabled(i));
		prog_name_suffix.push_back(enabled ? '1' : '0');
		name[12] = char('0'+i);
		set_bool_prefixes(name, enabled, shaders_enabled);
	}
}


void shader_t::setup_scene_bounds() const {

	float const scene_zmin(get_zval_min()), scene_zmax(get_zval_max());
	add_uniform_vector3d("scene_llc",   vector3d(-X_SCENE_SIZE, -Y_SCENE_SIZE, scene_zmin));
	add_uniform_vector3d("scene_scale", vector3d(2.0*X_SCENE_SIZE, 2.0*Y_SCENE_SIZE, (scene_zmax - scene_zmin)));
}


void shader_t::setup_fog_scale() const {

	add_uniform_float("fog_scale", (fog_enabled ? 1.0 : 0.0));
}


void shader_t::check_for_fog_disabled() {

	if (!fog_enabled) {for (unsigned d = 0; d < 2; ++d) {set_prefix("#define NO_FOG", d);}} // VS/FS
}


void shader_t::set_prefix(char const *const prefix, unsigned shader_type) {

	assert(shader_type < NUM_SHADER_TYPES);
	prog_name_suffix.push_back(',');
	prog_name_suffix.push_back('s');
	prog_name_suffix.push_back('0'+shader_type);
	prog_name_suffix += prefix;
	prepend_string[shader_type] += prefix;
	prepend_string[shader_type].push_back('\n');
}


void shader_t::set_bool_prefix(char const *const name, bool val, unsigned shader_type) {
	
	set_prefix_str((string("const bool ") + name + (val ? " = true;" : " = false;")), shader_type);
}


void shader_t::set_bool_prefixes(char const *const name, bool val, unsigned shaders_enabled) {

	string const prefix(string("const bool ") + name + (val ? " = true;" : " = false;"));

	for (unsigned s = 0; s < NUM_SHADER_TYPES; ++s) { // put into correct shader(s): V, F, G, TC, TE
		if (shaders_enabled & (1<<s)) {set_prefix_str(prefix, s);}
	}
}


void shader_t::set_int_prefix(char const *const name, int val, unsigned shader_type) {
	
	ostringstream oss;
	oss << val;
	set_prefix_str((string("const int ") + name + " = " + oss.str() + ";"), shader_type);
}


// *** shader and program setup ***


struct program_t {
	unsigned p, sixs[NUM_SHADER_TYPES];

	program_t() : p(0) {
		for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {sixs[i] = 0;}
	}
	program_t(unsigned p_, unsigned sixs_[NUM_SHADER_TYPES]) : p(p_) {
		for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {sixs[i] = sixs_[i];}
	}
};


class string_prog_map : public map<string, program_t> {
public:
	void clear() {
		free_data();
		map<string, program_t>::clear();
	}
	void free_data() {
		for (const_iterator s = begin(); s != end(); ++s) {
			for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
				if (s->second.sixs[i]) {glDetachShader(s->second.p, s->second.sixs[i]);}
			}
			glDeleteProgram(s->second.p);
		}
	}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~string_prog_map() {free_data();}
};


class string_shad_map : public map<string, unsigned> {
public:
	void clear() {
		free_data();
		map<string, unsigned>::clear();
	}
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			//assert(glIsShader(i->second));
			glDeleteShader(i->second);
		}
	}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~string_shad_map() {free_data();}
};


string_prog_map loaded_programs;
string_shad_map loaded_shaders[NUM_SHADER_TYPES]; // vertex=0, fragment=1, geometry=2, tess_control=3, tess_eval=4
map<string, string> loaded_files;


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
	cout << "GLSL Shader Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader && GL_EXT_geometry_shader4) return 1;
	cerr << "Error setting up vertex and fragment GLSL shaders." << endl;
	return 0;
}


void clear_shaders() {

	loaded_programs.clear();

	for (unsigned d = 0; d < NUM_SHADER_TYPES; ++d) {
		loaded_shaders[d].clear();
	}
	clear_cached_shaders();
}


void reload_all_shaders() { // clears and reloads *everything*

	// Note: do we want/need some function called every frame that check if shader files have been modified and calls this?
	cout << "Reloading all shaders" << endl;
	clear_shaders();
	loaded_files.clear();
}


bool shader_t::load_shader_file(string const &fname, string &data) {

	if (fname.empty()) return 0;
	map<string, string>::const_iterator i(loaded_files.find(fname));
	
	if (i != loaded_files.end()) {
		data += i->second;
		return 1;
	}
	ifstream in(fname.c_str());
	if (!in.good()) return 0;
	string line, file_contents;
	while (std::getline(in, line)) {file_contents += line + '\n';}
	loaded_files[fname] = file_contents;
	if (PRINT_SHADER) cout << "shader data:" << endl << file_contents << endl;
	data += file_contents;
	return 1;
}


bool shader_t::clear_shader_file(string const &fname) {

	assert(!fname.empty());
	map<string, string>::iterator i(loaded_files.find(fname));
	if (i == loaded_files.end()) return 0; // not found - error?
	loaded_files.erase(i);
	return 1;
}


void shader_t::filename_split(string const &fname, vector<string> &fns, char sep) {

	stringstream ss(fname);
    string fn;
	while (getline(ss, fn, sep)) {fns.push_back(fn);}
}


void shader_t::get_shader_filenames(string const &name, unsigned type, vector<string> &fns) {

	string const shader_name_table[NUM_SHADER_TYPES] = {"vert", "frag", "geom", "tess_control", "tess_eval"};
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


unsigned shader_t::get_shader(string const &name, unsigned type) const {
	
	//RESET_TIME;
	if (name.empty()) return 0; // none selected
	int const shader_type_table[NUM_SHADER_TYPES] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_GEOMETRY_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER};
	assert(type < NUM_SHADER_TYPES);
	string const lookup_name(name + prepend_string[type]);
	string_shad_map::const_iterator it(loaded_shaders[type].find(lookup_name));
	if (it != loaded_shaders[type].end()) {return it->second;} // already loaded
	
	// create a new shader
	string const version_info("#version 400\n");
	vector<string> fns;
	get_shader_filenames(name, type, fns);
	bool failed(0);
	unsigned shader(0);

	while (1) { // retry loop
		if (failed) {
			if (!yes_no_query("Retry?")) {
				cerr << "Exiting." << endl;
				exit(1);
			}
			for (vector<string>::const_iterator i = fns.begin(); i != fns.end(); ++i) {
				if (clear_shader_file(*i)) {cout << "Reloading shader component " << *i << endl;}
			}
			failed = 0;
		}
		string data(version_info + prepend_string[type]);

		for (vector<string>::const_iterator i = fns.begin(); i != fns.end(); ++i) {
			if (!load_shader_file(*i, data)) {
				cerr << "Error loading shader file " << *i << "." << endl;
				failed = 1; break;
			}
		}
		if (failed) continue;
		if (DEBUG_SHADER) {cout << "final shader data for <" << name << ">:" << endl << data << endl;}
		shader = glCreateShader(shader_type_table[type]);
	
		if (!shader) {
			cerr << "Error: Failed to create shader " << name << "." << endl;
			failed = 1; continue;
		}
		const char *src(data.c_str());
		glShaderSource(shader, 1, &src, 0);
		glCompileShader(shader);
		int status(0);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

		if (status != GL_TRUE) {
			cerr << "Compilation of shader " << name << " failed with status " << status << endl;
			print_shader_info_log(shader);
			cerr << endl;
			failed = 1; continue;
		}
		break;
	} // while(1)
	if (PRINT_LOG) {print_shader_info_log(shader);}
	loaded_shaders[type][lookup_name] = shader; // cache the shader
	//PRINT_TIME("Create Shader");
	return shader;
}


// See http://www.lighthouse3d.com/tutorials/glsl-core-tutorial/glsl-core-tutorial-create-a-program/
bool shader_t::begin_shader(bool do_enable) {

	//RESET_TIME;
	// get the program
	string pname(prog_name_suffix);
	for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {pname += shader_names[i] + ",";} // unique program identifier
	string_prog_map::const_iterator it(loaded_programs.find(pname));
	program = 0;

	if (it != loaded_programs.end()) { // program already exists
		program = it->second.p;
	}
	else { // create a new program
		program = glCreateProgram();
		unsigned shader_ixs[NUM_SHADER_TYPES];
		
		for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
			// Note: we *can* attach multiple shaders of the same type to a single program as long as only one has a main()
			shader_ixs[i] = get_shader(shader_names[i], i);
			if (shader_ixs[i]) {glAttachShader(program, shader_ixs[i]);}
		}
		assert(shader_ixs[0] && shader_ixs[1]); // vertex and fragment shaders are required, geometry shader is optional

		if (shader_ixs[2]) { // setup geometry shader
			assert(GL_EXT_geometry_shader4);
			glProgramParameteriEXT(program, GL_GEOMETRY_INPUT_TYPE_EXT,  in_prim);
			glProgramParameteriEXT(program, GL_GEOMETRY_OUTPUT_TYPE_EXT, out_prim);
			int max_verts_out(0);
			glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &max_verts_out); // get max
			if (verts_out == 0) {verts_out = max_verts_out;} else {assert(verts_out <= max_verts_out);}
			assert(verts_out > 0);
			glProgramParameteriEXT(program, GL_GEOMETRY_VERTICES_OUT_EXT, verts_out);
		}
		check_gl_error(300);
		glLinkProgram(program);
		int status(0);
		glGetProgramiv(program, GL_LINK_STATUS, &status);

		if (status != GL_TRUE) {
			cerr << "Linking of program " << pname << " failed with status " << status << endl;
			print_program_info_log();
			cout << endl;
			exit(1);
		}
		if (PRINT_LOG) print_program_info_log();
		loaded_programs[pname] = program_t(program, shader_ixs); // cache the program
		//PRINT_TIME("Create Program"); // 90ms
	}
	if (do_enable) {enable();}
	return 1;
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
	program = 0;
	
	for (unsigned i = 0; i < NUM_SHADER_TYPES; ++i) {
		prepend_string[i].clear();
		shader_names[i].clear();
	}
	prog_name_suffix.clear();
	attrib_locs.clear();
}


// some simple shared shaders
void shader_t::begin_color_only_shader() {

	set_vert_shader("vert_xform_only");
	set_frag_shader("color_only");
	begin_shader();
}

void shader_t::begin_simple_textured_shader(float min_alpha, bool include_2_lights, bool use_texgen) {

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
	add_uniform_float("min_alpha", min_alpha);
	add_uniform_int("tex0", 0);
	if (use_fog) {add_uniform_float("fog_scale", 1.0);}
}

void shader_t::begin_untextured_lit_glcolor_shader() {

	set_prefix("#define USE_LIGHT_COLORS", 0); // VS
	begin_simple_textured_shader(0.0, 1); // lighting (not actually textured)
	select_texture(WHITE_TEX); // untextured
}


// **************** INSTANCING ****************


void instance_render_t::add_cur_inst() {

	xform_matrix xf;
	xf.assign_mv_from_gl();
	add_inst(xf);
}


// Note: assumes cur_vbo is currently bound by the caller, and will leave it bound after the call
void instance_render_t::draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo, int index_type, void *indices) { // indices can be NULL

	if (inst_xforms.empty()) return;

	if (loc < 0) { // hardware instancing not used, so we iterate
		glPushMatrix();

		for (vector<xform_matrix>::const_iterator i = inst_xforms.begin(); i != inst_xforms.end(); ++i) {
			glLoadIdentity();
			i->apply();
			
			if (index_type != GL_NONE) { // indexed
				glDrawElements(prim_type, count, index_type, indices);
			}
			else {
				assert(indices == NULL);
				glDrawArrays(prim_type, 0, count); // hard-coded first=0
			}
		}
		glPopMatrix();
	}
	else { // use hardware instancing
		// we need to upload the transforms here, but if there is a current vbo associated with the object data
		// we need to unbind it, upload the transforms, then rebind the original vbo
		if (cur_vbo) {bind_vbo(0);}
		shader_float_matrix_uploader<4>::enable(loc, 1, inst_xforms.front().get_ptr());
		if (cur_vbo) {bind_vbo(cur_vbo);}
	
		if (index_type != GL_NONE) { // indexed
			glDrawElementsInstanced(prim_type, count, index_type, indices, inst_xforms.size());
		}
		else {
			assert(indices == NULL);
			glDrawArraysInstanced(prim_type, 0, count, inst_xforms.size()); // hard-coded first=0
		}
		shader_float_matrix_uploader<4>::disable(loc);
	}
	inst_xforms.clear();
}

