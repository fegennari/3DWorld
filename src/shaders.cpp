// 3D World - Vertex and Fragment GLSL Shader Framework
// by Frank Gennari
// 11/1/10
#include "GL/glew.h"
#include "3DWorld.h"
#include <fstream>

using namespace std;


bool const PRINT_SHADER = 0;
bool const PRINT_LOG    = 0;

string const shaders_dir = "shaders";


// *** uniform variables setup ***


typedef pair<float const *, unsigned> float_arr_t;
typedef map<string, float_arr_t> u_float_array_map_t;
typedef map<string, float> u_float_map_t;
typedef map<string, int  > u_int_map_t;
u_float_array_map_t u_float_array_map;
u_float_map_t u_float_map;
u_int_map_t   u_int_map;
string prepend_string[3]; // vertex=0, fragment=1, geometry=3
string prog_name_suffix;


void add_uniform_float_array(string const &name, float const *const val, unsigned num) {

	assert(!name.empty());
	u_float_array_map[name] = make_pair(val, num);
}


void add_uniform_float(string const &name, float val) {

	assert(!name.empty());
	u_float_map[name] = val;
}


void add_uniform_int(string const &name, int val) {

	assert(!name.empty());
	u_int_map[name] = val;
}


int get_uniform_loc(int program, string const &name) {

	int const loc(glGetUniformLocation(program, name.c_str()));
	//cout << "name: " << name << ", loc: " << loc << endl;
	//assert(loc >= 0); // Note: if variable is unused, loc will be -1
	return loc;
}


bool set_uniform_buffer(int program) {

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

	//SurfaceColor might change, so we'll query its offset/size.
	const char *name = "data";
	unsigned index;
	int offset, size;
	//First, get the index for the uniform
	glGetUniformIndices(program, 1, &name, &index);
	assert((int)index >= 0);
	//Use the index to query offset and size
	glGetActiveUniformsiv(program, 1, &index, GL_UNIFORM_OFFSET, &offset);
	glGetActiveUniformsiv(program, 1, &index, GL_UNIFORM_SIZE, &size);
	//cout << "index: " << index << ", offset: " << offset << ", size: " << size << endl;

	// Create UBO
	static unsigned buffer_id(0);
	if (buffer_id == 0 || !glIsBuffer(buffer_id)) glGenBuffers(1, &buffer_id);
	assert(buffer_id > 0);
	glBindBuffer(GL_UNIFORM_BUFFER, buffer_id);
	
	// We can use BufferData to upload our data to the shader, since we know it's in the std140 layout
	glBufferData(GL_UNIFORM_BUFFER, uniformBlockSize, NULL, GL_DYNAMIC_DRAW);

	// Bind constants to UBO binding point 0
	char data[1000] = {0};
	glBindBufferBase(GL_UNIFORM_BUFFER, 0, buffer_id);
	glBufferSubData(GL_UNIFORM_BUFFER, offset, size, &data);
	glBindBuffer(GL_UNIFORM_BUFFER, 0); // unbind
	return 1;
}


void setup_uniforms(int program) {

	assert(program);
	set_uniform_buffer(program);

	for (u_float_array_map_t::const_iterator i = u_float_array_map.begin(); i != u_float_array_map.end(); ++i) {
		int const loc(get_uniform_loc(program, i->first));
		if (loc >= 0) glUniform1fv(loc, i->second.second, i->second.first);
	}
	for (u_float_map_t::const_iterator i = u_float_map.begin(); i != u_float_map.end(); ++i) {
		int const loc(get_uniform_loc(program, i->first));
		if (loc >= 0) glUniform1f(loc, i->second);
	}
	for (u_int_map_t::const_iterator i = u_int_map.begin(); i != u_int_map.end(); ++i) {
		int const loc(get_uniform_loc(program, i->first));
		if (loc >= 0) glUniform1i(loc, i->second);
	}
}


void setup_enabled_lights(unsigned num) {

	prog_name_suffix += ",el";

	for (unsigned i = 0; i < num; ++i) { // 0=sun, 1=moon
		GLboolean const enabled(glIsEnabled(GL_LIGHT0 + i));
		prog_name_suffix += (enabled ? '1' : '0');

		for (unsigned s = 0; s < 2; ++s) { // put into vertex and fragment shaders
			set_bool_shader_prefix((string("enable_light") + char('0'+i)), (enabled != 0), s);
		}
	}
}


void setup_fog_scale() {

	add_uniform_float("fog_scale", (glIsEnabled(GL_FOG) ? 1.0 : 0.0));
}


void set_shader_prefix(string const &prefix, unsigned shader_type) {

	assert(shader_type < 3);
	prog_name_suffix += ",s" + ('0'+shader_type) + prefix;
	prepend_string[shader_type] += prefix + '\n';
}


void set_bool_shader_prefix(string const &name, bool val, unsigned shader_type) {
	
	set_shader_prefix((string("const bool ") + name + " = " + (val ? "true;" : "false;")), shader_type);
}


// *** shader and program setup ***


struct program_t {
	unsigned p, vs, fs, gs;
	program_t(unsigned p_=0, unsigned vs_=0, unsigned fs_=0, unsigned gs_=0) : p(p_), vs(vs_), fs(fs_), gs(gs_) {}
};


class string_prog_map : public map<string, program_t> {
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			if (i->second.vs) glDetachShader(i->second.p, i->second.vs);
			if (i->second.fs) glDetachShader(i->second.p, i->second.fs);
			if (i->second.gs) glDetachShader(i->second.p, i->second.gs);
			glDeleteProgram(i->second.p);
		}
	}
public:
	~string_prog_map() {free_data();}
};


class string_shad_map : public map<string, unsigned> {
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			glDeleteShader(i->second);
		}
	}
public:
	~string_shad_map() {free_data();}
};


string_prog_map loaded_programs;
string_shad_map loaded_shaders[3]; // vertex=0, fragment=1, geometry=3
map<string, string> loaded_files;


bool load_shader_file(string const &fname, string &data) {

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


bool setup_shaders() {

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


void filename_split(string const &fname, vector<string> &fns, char sep) {

	stringstream ss(fname);
    string fn;
    while(getline(ss, fn, sep)) fns.push_back(fn);
}


void print_shader_info_log(unsigned shader) {

	int len(0), len2(0);
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);

	if (len > 0) {
		vector<char> info_log_msg(len);
		glGetShaderInfoLog(shader, len, &len2, &info_log_msg.front()); 
		assert(len2 <= len);
		cout << "Info log: " << string(info_log_msg.begin(), info_log_msg.end()) << endl;
	}
}


unsigned get_shader(string const &name, unsigned type) {
	
	RESET_TIME;
	int const shader_type_table   [3] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_GEOMETRY_SHADER_EXT};
	string const shader_name_table[3] = {"vert", "frag", "geom"};
	assert(type < 3);
	if (name.empty()) return 0; // none selected
	string const lookup_name(name + prepend_string[type]);
	string_shad_map::const_iterator it(loaded_shaders[type].find(lookup_name));
	if (it != loaded_shaders[type].end()) return it->second; // already loaded
	
	// create a new shader
	string data(prepend_string[type]);
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
		exit(1);
	}
	if (PRINT_LOG) print_shader_info_log(shader);
	loaded_shaders[type][lookup_name] = shader; // cache the shader
	//PRINT_TIME("Create Shader"); // 43ms
	return shader;
}


void print_program_info_log(unsigned program) {

	int len(0), len2(0);
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);

	if (len > 0) {
		vector<char> info_log_msg(len);
		glGetProgramInfoLog(program, len, &len2, &info_log_msg.front()); 
		assert(len2 <= len);
		cout << "Info log: " << string(info_log_msg.begin(), info_log_msg.end()) << endl;
	}
}


bool set_shader_prog(string const &vs_name, string const &fs_name, string const &gs_name,
					 int in_prim, int out_prim, int verts_out)
{
	// get the program
	RESET_TIME;
	string const pname(vs_name + "," + fs_name + "," + gs_name + "," + prog_name_suffix); // unique program identifier
	string_prog_map::const_iterator it(loaded_programs.find(pname));
	unsigned program(0);

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
			if (verts_out == 0) glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &verts_out); // get max
			assert(verts_out > 0);
			glProgramParameteriEXT(program, GL_GEOMETRY_VERTICES_OUT_EXT, verts_out);
		}
		glLinkProgram(program);
		int status(0);
		glGetProgramiv(program, GL_LINK_STATUS, &status);

		if (status != GL_TRUE) {
			cerr << "Linking of program " << pname << " failed with status " << status << endl;
			print_program_info_log(program);
			exit(1);
		}
		if (PRINT_LOG) print_program_info_log(program);
		loaded_programs[pname] = program_t(program, vs, fs, gs); // cache the program
		//PRINT_TIME("Create Program"); // 90ms
	}
	assert(program);
	glUseProgram(program);
	setup_uniforms(program);
	return 1; // can't fail yet
}


void unset_shader_prog() {

	glUseProgram(0);
	u_float_array_map.clear();
	u_float_map.clear();
	u_int_map.clear();
	
	for (unsigned i = 0; i < 3; ++i) {
		prepend_string[i].clear();
	}
	prog_name_suffix.clear();
}

