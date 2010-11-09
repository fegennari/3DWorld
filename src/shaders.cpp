// 3D World - Vertex and Fragment GLSL Shader Framework
// by Frank Gennari
// 11/1/10
#include "GL/glew.h"
#include "3DWorld.h"
#include <fstream>

using namespace std;


string const shaders_dir = "shaders";


// uniform variables setup


void setup_enabled_lights(int program) {

	int enabled_gl_lights(0);

	for (unsigned i = 0; i < 8; ++i) { // max of 8 lights (GL_LIGHT0 - GL_LIGHT7): sun, moon, lightning
		int const light(GL_LIGHT0 + i); // should be sequential
		if (glIsEnabled(light)) enabled_gl_lights |= (1 << i);
	}
	int const loc1(glGetUniformLocation(program, "enabled_gl_lights"));
	glUniform1i(loc1, enabled_gl_lights);
}


struct program_t {
	unsigned p, vs, fs;
	program_t(unsigned p_=0, unsigned vs_=0, unsigned fs_=0) : p(p_), vs(vs_), fs(fs_) {}
};


class string_prog_map : public map<string, program_t> {
	void free_data() {
		for (const_iterator i = begin(); i != end(); ++i) {
			glDetachShader(i->second.p, i->second.vs);
			glDetachShader(i->second.p, i->second.fs);
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
string_shad_map loaded_shaders[2]; // vertex=0, fragment=1


bool load_shader_file(string const &fname, string &data) {

	data.clear();
	if (fname.empty()) return 0;
	ifstream in(fname.c_str());
	if (!in.good()) return 0;
	string line;
	while (std::getline(in, line)) data += line + '\n';
	//cout << "shader data:" << endl << data << endl;
	return 1;
}


bool setup_shaders() {

	if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader) return 1;
	cerr << "Error setting up vertex and fragment GLSL shaders." << endl;
	return 0;
}


void clear_shaders() {

	loaded_programs.clear();

	for (unsigned d = 0; d < 2; ++d) {
		loaded_shaders[d].clear();
	}
}


unsigned get_shader(string const &name, bool vf) {
	
	assert(!name.empty());
	string_shad_map::const_iterator it(loaded_shaders[vf].find(name));
	if (it != loaded_shaders[vf].end()) return it->second; // already loaded

	// create a new shader
	string const fname(shaders_dir + "/" + name + (vf ? ".frag" : ".vert"));
	string data;
	
	if (!load_shader_file(fname, data)) {
		cerr << "Error loading shader file " << fname << ". Exiting." << endl;
		exit(1);
	}
	unsigned const shader(glCreateShader(vf ? GL_FRAGMENT_SHADER : GL_VERTEX_SHADER));
	assert(shader);
	const char *src(data.c_str());
	glShaderSource(shader, 1, &src, 0);
	glCompileShader(shader);
	int status(0);
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

	if (status != GL_TRUE) {
		cerr << "Compilation of shader " << name << " failed with status " << status << endl;
		int len(0), len2(0);
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);

		if (len > 0) {
			vector<char> info_log_msg(len);
			glGetShaderInfoLog(shader, len, &len2, &info_log_msg.front()); 
			assert(len2 <= len);
			cerr << "Info log: " << string(info_log_msg.begin(), info_log_msg.end()) << endl;
		}
		exit(1);
	}
	loaded_shaders[vf][name] = shader; // cache the shader
	return shader;
}


bool set_shader_prog(string const &vs_name, string const &fs_name) {

	// get the program
	string const pname(vs_name + "," + fs_name); // unique program identifier
	string_prog_map::const_iterator it(loaded_programs.find(pname));
	unsigned program(0);

	if (it != loaded_programs.end()) { // program already exists
		program = it->second.p;
	}
	else { // create a new program
		program = glCreateProgram();
		unsigned const vs(get_shader(vs_name, 0));
		unsigned const fs(get_shader(fs_name, 1));
		glAttachShader(program, vs);
		glAttachShader(program, fs);
		glLinkProgram(program);
		int status(0);
		glGetProgramiv(program, GL_LINK_STATUS, &status);

		if (status != GL_TRUE) {
			cerr << "Linking of program " << pname << " failed with status " << status << endl;
			int len(0), len2(0);
			glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);

			if (len > 0) {
				vector<char> info_log_msg(len);
				glGetProgramInfoLog(program, len, &len2, &info_log_msg.front()); 
				assert(len2 <= len);
				cerr << "Info log: " << string(info_log_msg.begin(), info_log_msg.end()) << endl;
			}
			exit(1);
		}
		loaded_programs[pname] = program_t(program, vs, fs); // cache the program
	}
	assert(program);
	glUseProgram(program);
	setup_enabled_lights(program);
	return 1; // can't fail yet
}


void unset_shader_prog() {

	glUseProgram(0);
}

