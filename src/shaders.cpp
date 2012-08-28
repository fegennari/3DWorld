// 3D World - Vertex and Fragment GLSL Shader Framework
// by Frank Gennari
// 11/1/10
#include "shaders.h"
#include "mesh.h" // for scene bounds
#include <fstream>

using namespace std;

bool const PRINT_SHADER = 0; // print shaders loaded from files
bool const DEBUG_SHADER = 0; // print final generated shaders
bool const PRINT_LOG    = 0;

string const shaders_dir = "shaders";

extern bool disable_shaders;


// *** uniform variables setup ***


char const *append_array_ix(string &s, unsigned i) {

	assert(i <= 9);
	s.push_back('[');
	s.push_back('0'+i);
	s.push_back(']');
	return s.c_str();
}


int shader_t::get_uniform_loc(char const *const name) const {

	if (disable_shaders) return 0;
	assert(program && name);
	int const loc(glGetUniformLocation(program, name));
	//cout << "name: " << name << ", loc: " << loc << endl;
	//assert(loc >= 0); // Note: if variable is unused, loc will be -1
	return loc;
}


void shader_t::set_uniform_float_array(int loc, float const *const val, unsigned num) const {
	if (loc >= 0) glUniform1fv(loc, num, val);
}

void shader_t::set_uniform_float(int loc, float val) const {
	if (loc >= 0) glUniform1f(loc, val);
}

void shader_t::set_uniform_int(int loc, int val) const {
	if (loc >= 0) glUniform1i(loc, val);
}

void shader_t::set_uniform_vector2d(int loc, vector2d const &val) const {
	if (loc >= 0) glUniform2fv(loc, 1, &val.x);
}

void shader_t::set_uniform_vector3d(int loc, vector3d const &val) const {
	if (loc >= 0) glUniform3fv(loc, 1, &val.x);
}

void shader_t::set_uniform_color(int loc, colorRGBA const &val) const {
	if (loc >= 0) glUniform4fv(loc, 1, &val.R);
}

void shader_t::set_uniform_color(int loc, colorRGB  const &val) const {
	if (loc >= 0) glUniform3fv(loc, 1, &val.R);
}

void shader_t::set_uniform_matrid_4x4(int loc, float *m, bool transpose) const {
	if (loc >= 0) glUniformMatrix4fv(loc, 1, transpose, m);
}


void shader_t::add_uniform_float_array(char const *const name, float const *const val, unsigned num) const {
	if (!disable_shaders) set_uniform_float_array(get_uniform_loc(name), val, num);
}

void shader_t::add_uniform_float(char const *const name, float val) const {
	if (!disable_shaders) set_uniform_float(get_uniform_loc(name), val);
}

void shader_t::add_uniform_int(char const *const name, int val) const {
	if (!disable_shaders) set_uniform_int(get_uniform_loc(name), val);
}

void shader_t::add_uniform_vector2d(char const *const name, vector2d const &val) const {
	if (!disable_shaders) set_uniform_vector2d(get_uniform_loc(name), val);
}

void shader_t::add_uniform_vector3d(char const *const name, vector3d const &val) const {
	if (!disable_shaders) set_uniform_vector3d(get_uniform_loc(name), val);
}

void shader_t::add_uniform_color(char const *const name, colorRGBA const &val) const {
	if (!disable_shaders) set_uniform_color(get_uniform_loc(name), val);
}

void shader_t::add_uniform_color(char const *const name, colorRGB  const &val) const {
	if (!disable_shaders) set_uniform_color(get_uniform_loc(name), val);
}

void shader_t::add_uniform_matrid_4x4(char const *const name, float *m, bool transpose) const {
	if (!disable_shaders) set_uniform_matrid_4x4(get_uniform_loc(name), m, transpose);
}


// unused, unfinished
bool shader_t::set_uniform_buffer_data(char const *name, float const *data, unsigned size, unsigned &buffer_id) const {

	if (disable_shaders) return 0;
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


void shader_t::register_attrib_name(char const *name, unsigned bind_ix) {

	if (disable_shaders) return;
	assert(bind_ix < 100); // sanity check
	int const loc(get_attrib_loc(name));
	if (bind_ix >= attrib_locs.size()) attrib_locs.resize(bind_ix+1);
	attrib_locs[bind_ix] = loc;
}


int shader_t::attrib_loc_by_ix(unsigned ix) const {

	assert(ix < attrib_locs.size());
	return attrib_locs[ix]; // is it legal for this to return -1?
}


void shader_t::add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const {

	if (disable_shaders) return;
	int const loc(attrib_loc_by_ix(ix));
	if (loc < 0) return;

	switch (num) {
		case 1: glVertexAttrib1fv(loc, val); break;
		case 2: glVertexAttrib2fv(loc, val); break;
		case 3: glVertexAttrib3fv(loc, val); break;
		case 4: glVertexAttrib4fv(loc, val); break;
		default: assert(0);
	}
}


void shader_t::add_attrib_float(unsigned ix, float val) const {

	if (disable_shaders) return;
	int const loc(attrib_loc_by_ix(ix));
	if (loc >= 0) glVertexAttrib1f(loc, val);
}


void shader_t::add_attrib_int(unsigned ix, int val) const {

	if (disable_shaders) return;
	int const loc(attrib_loc_by_ix(ix));
	if (loc >= 0) glVertexAttrib1s(loc, val);
}


// *** other variables setup ***


void shader_t::setup_enabled_lights(unsigned num, unsigned shaders_enabled) {

	assert(num <= 8);
	prog_name_suffix += ",el";
	string name("enable_light0");

	for (unsigned i = 0; i < num; ++i) { // 0=sun, 1=moon, ...
		GLboolean const enabled(glIsEnabled(GL_LIGHT0 + i));
		prog_name_suffix += (enabled ? '1' : '0');
		name.back() = char('0'+i);

		for (unsigned s = 0; s < 3; ++s) { // put into correct shader(s): V, F, G
			if (shaders_enabled & (1<<s)) {set_bool_prefix(name, (enabled != 0), s);}
		}
	}
}


void shader_t::setup_scene_bounds() const {

	float const scene_zmin(get_zval_min()), scene_zmax(get_zval_max());
	add_uniform_vector3d("scene_llc",   vector3d(-X_SCENE_SIZE, -Y_SCENE_SIZE, scene_zmin));
	add_uniform_vector3d("scene_scale", vector3d(2.0*X_SCENE_SIZE, 2.0*Y_SCENE_SIZE, (scene_zmax - scene_zmin)));
}


void shader_t::setup_fog_scale() const {

	add_uniform_float("fog_scale", (glIsEnabled(GL_FOG) ? 1.0 : 0.0));
}


void shader_t::set_prefix(string const &prefix, unsigned shader_type) {

	assert(shader_type < 3);
	prog_name_suffix += ",s" + ('0'+shader_type) + prefix;
	prepend_string[shader_type] += prefix + '\n';
}


void shader_t::set_bool_prefix(string const &name, bool val, unsigned shader_type) {
	
	set_prefix((string("const bool ") + name + " = " + (val ? "true;" : "false;")), shader_type);
}


void shader_t::set_int_prefix(string const &name, int val, unsigned shader_type) {
	
	ostringstream oss;
	oss << val;
	set_prefix((string("const int ") + name + " = " + oss.str() + ";"), shader_type);
}


// *** shader and program setup ***


struct program_t {
	unsigned p, vs, fs, gs;
	program_t(unsigned p_=0, unsigned vs_=0, unsigned fs_=0, unsigned gs_=0) : p(p_), vs(vs_), fs(fs_), gs(gs_) {}
};


class string_prog_map : public map<string, program_t> {
public:
	void clear() {
		free_data();
		map<string, program_t>::clear();
	}
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			if (i->second.vs) glDetachShader(i->second.p, i->second.vs);
			if (i->second.fs) glDetachShader(i->second.p, i->second.fs);
			if (i->second.gs) glDetachShader(i->second.p, i->second.gs);
			glDeleteProgram(i->second.p);
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
string_shad_map loaded_shaders[3]; // vertex=0, fragment=1, geometry=3
map<string, string> loaded_files;


bool setup_shaders() {

	cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL Shader Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader && GL_EXT_geometry_shader4) return 1;
	cerr << "Error setting up vertex and fragment GLSL shaders." << endl;
	return 0;
}


void clear_shaders() {

	loaded_programs.clear();

	for (unsigned d = 0; d < 3; ++d) {
		loaded_shaders[d].clear();
	}
}


bool shader_t::load_shader_file(string const &fname, string &data) const {

	if (fname.empty()) return 0;
	map<string, string>::const_iterator i(loaded_files.find(fname));
	
	if (i != loaded_files.end()) {
		data += i->second;
		return 1;
	}
	ifstream in(fname.c_str());
	if (!in.good()) return 0;
	string line, file_contents;
	while (std::getline(in, line)) file_contents += line + '\n';
	loaded_files[fname] = file_contents;
	if (PRINT_SHADER) cout << "shader data:" << endl << file_contents << endl;
	data += file_contents;
	return 1;
}


void shader_t::filename_split(string const &fname, vector<string> &fns, char sep) const {

	stringstream ss(fname);
    string fn;
    while (getline(ss, fn, sep)) fns.push_back(fn);
}


unsigned shader_t::get_shader(string const &name, unsigned type) const {
	
	//RESET_TIME;
	if (name.empty()) return 0; // none selected
	int const shader_type_table   [3] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_GEOMETRY_SHADER_EXT};
	string const shader_name_table[3] = {"vert", "frag", "geom"};
	assert(type < 3);
	string const lookup_name(name + prepend_string[type]);
	string_shad_map::const_iterator it(loaded_shaders[type].find(lookup_name));
	if (it != loaded_shaders[type].end()) return it->second; // already loaded
	
	// create a new shader
	string const version_info("#version 400\n");
	string data(version_info + prepend_string[type]);
	vector<string> fns;
	filename_split(name, fns, '+');

	for (vector<string>::const_iterator i = fns.begin(); i != fns.end(); ++i) {
		assert(!i->empty());
		string fname(shaders_dir + "/" + *i);
		
		if ((*i)[i->size()-1] == '*') { // wildcard shader file: works with all shader types
			fname.erase(fname.size()-1);
		}
		else { // add shader type extension
			fname += "." + shader_name_table[type];
		}
		if (!load_shader_file(fname, data)) {
			cerr << "Error loading shader file " << fname << ". Exiting." << endl;
			exit(1);
		}
	}
	if (DEBUG_SHADER) cout << "final shader data for <" << name << ">:" << endl << data << endl;
	unsigned const shader(glCreateShader(shader_type_table[type]));
	assert(shader);
	const char *src(data.c_str());
	glShaderSource(shader, 1, &src, 0);
	glCompileShader(shader);
	int status(0);
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

	if (status != GL_TRUE) {
		cerr << "Compilation of shader " << name << " failed with status " << status << endl;
		print_shader_info_log(shader);
		cout << endl;
		exit(1);
	}
	if (PRINT_LOG) print_shader_info_log(shader);
	loaded_shaders[type][lookup_name] = shader; // cache the shader
	//PRINT_TIME("Create Shader"); // 43ms
	return shader;
}


bool shader_t::begin_shader() {

	if (disable_shaders) return 0;
	// get the program
	//RESET_TIME;
	string const pname(vs_name + "," + fs_name + "," + gs_name + "," + prog_name_suffix); // unique program identifier
	string_prog_map::const_iterator it(loaded_programs.find(pname));
	program = 0;

	if (it != loaded_programs.end()) { // program already exists
		program = it->second.p;
	}
	else { // create a new program
		program = glCreateProgram();
		unsigned const vs(get_shader(vs_name, 0));
		unsigned const fs(get_shader(fs_name, 1));
		unsigned const gs(get_shader(gs_name, 2));
		assert(vs && fs); // vertex and fragment shaders are required, geometry shader is optional
		if (vs) glAttachShader(program, vs);
		if (fs) glAttachShader(program, fs);
		if (gs) glAttachShader(program, gs);

		if (gs) { // setup geometry shader
			// Note: we MUST *NOT* be in a display list when we get here
			assert(GL_EXT_geometry_shader4);
			glProgramParameteriEXT(program, GL_GEOMETRY_INPUT_TYPE_EXT, in_prim);
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
		loaded_programs[pname] = program_t(program, vs, fs, gs); // cache the program
		//PRINT_TIME("Create Program"); // 90ms
	}
	assert(program);
	glUseProgram(program);
	return 1;
}


void shader_t::print_shader_info_log(unsigned shader) const {

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

	glUseProgram(0);
	program = 0;
	
	for (unsigned i = 0; i < 3; ++i) {
		prepend_string[i].clear();
	}
	prog_name_suffix.clear();
	attrib_locs.clear();
	vs_name.clear();
	fs_name.clear();
	gs_name.clear();
}

