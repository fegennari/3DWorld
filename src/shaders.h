// 3D World - Shader Class Declaration
// by Frank Gennari
// 8/4/10

#ifndef _SHADERS_H_
#define _SHADERS_H_

#include "3DWorld.h"
#include "transform_obj.h"

using std::string;


unsigned const TEX0_S_ATTR = 0;
unsigned const TEX0_T_ATTR = 1;
unsigned const NUM_SHADER_TYPES = 5;


class shader_t {

	unsigned program; // active program
	string prepend_string[NUM_SHADER_TYPES]; // vertex=0, fragment=1, geometry=2, tess_control=3, tess_eval=4
	string prog_name_suffix;
	vector<int> attrib_locs;

	string shader_names[NUM_SHADER_TYPES];
	int in_prim, out_prim, verts_out;

	void filename_split(string const &fname, vector<string> &fns, char sep) const;
	bool load_shader_file(string const &fname, string &data) const;
	unsigned get_shader(string const &name, unsigned type) const;
	void print_shader_info_log(unsigned shader) const;
	void print_program_info_log() const;

public:
	shader_t() : program(0), in_prim(0), out_prim(0), verts_out(0) {}
	//~shader_t() {assert(!program);} // end_shader() should have been called (but not for cached global variables)

	void set_vert_shader(string const &vs_name_) {shader_names[0] = vs_name_;}
	void set_frag_shader(string const &fs_name_) {shader_names[1] = fs_name_;}
	void set_geom_shader(string const &gs_name_, int in_prim_, int out_prim_, int verts_out_) {
		shader_names[2] = gs_name_; in_prim = in_prim_; out_prim = out_prim_; verts_out = verts_out_;
	}
	void set_tess_control_shader(string const &tcs_name_) {shader_names[3] = tcs_name_;}
	void set_tess_eval_shader   (string const &tes_name_) {shader_names[4] = tes_name_;}
	void init_from_active_shader(shader_t const &s) {assert(s.is_setup()); program = s.program;}
	bool is_setup() const {return (program > 0);}
	void enable () const {assert(program); glUseProgram(program);}
	static void disable() {glUseProgram(0);}
	bool begin_shader(bool do_enable=1);
	void end_shader();
	void begin_color_only_shader();

	int get_uniform_loc(char const *const name) const;
	static void set_uniform_float_array(int loc, float const *const val, unsigned num);
	static void set_uniform_float      (int loc, float val);
	static void set_uniform_int        (int loc, int val);
	static void set_uniform_vector2d   (int loc, vector2d const &val);
	static void set_uniform_vector3d   (int loc, vector3d const &val);
	static void set_uniform_color      (int loc, colorRGBA const &val);
	static void set_uniform_color      (int loc, colorRGB  const &val);
	static void set_uniform_matrid_4x4 (int loc, float *m, bool transpose);

	void add_uniform_float_array (char const *const name, float const *const val, unsigned num) const;
	void add_uniform_float       (char const *const name, float val) const;
	void add_uniform_int         (char const *const name, int val) const;
	void add_uniform_vector2d    (char const *const name, vector2d const &val) const;
	void add_uniform_vector3d    (char const *const name, vector3d const &val) const;
	void add_uniform_color       (char const *const name, colorRGBA const &val) const;
	void add_uniform_color       (char const *const name, colorRGB  const &val) const;
	void add_uniform_matrid_4x4  (char const *const name, float *m, bool transpose) const;
	bool set_uniform_buffer_data (char const *name, float const *data, unsigned size, unsigned &buffer_id) const;

	int attrib_loc_by_ix(unsigned ix) const;
	int get_attrib_loc(char const *const name, bool allow_fail=0) const;
	void register_attrib_name(char const *const name, unsigned bind_ix);
	void set_attrib_float_array(int loc, float const *const val, unsigned num) const;
	void set_attrib_float      (int loc, float val) const;
	void set_attrib_int        (int loc, int val) const;
	void add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const;
	void add_attrib_float      (unsigned ix, float val) const;
	void add_attrib_int        (unsigned ix, int val) const;

	void setup_enabled_lights(unsigned num=2, unsigned shaders_enabled=3);
	void setup_scene_bounds() const;
	void setup_fog_scale() const;
	void check_for_fog_disabled();
	void set_prefix_str(string const &prefix, unsigned shader_type) {set_prefix(prefix.c_str(), shader_type);}
	void set_prefix(char const *const prefix, unsigned shader_type);
	void set_bool_prefix(char const *const name, bool val, unsigned shader_type);
	void set_bool_prefixes(char const *const name, bool val, unsigned shaders_enabled=3);
	void set_int_prefix(char const *const name, int val, unsigned shader_type);
};


template<unsigned N> struct shader_float_matrix_uploader {

	static void enable(int start_loc, int divisor, float const *const data=NULL) {
		assert(start_loc >= 0 && divisor >= 0);

		for (unsigned i = 0; i < N; ++i) {
			int const loc(start_loc + i);
			glEnableVertexAttribArray(loc);
			glVertexAttribPointer(loc, N, GL_FLOAT, GL_FALSE, N*N*sizeof(float), (const void *)(data + N*i));
			glVertexAttribDivisor(loc, divisor);
		}
	}
	static void disable(int start_loc) {
		for (unsigned i = 0; i < N; ++i) {
			int const loc(start_loc + i);
			glVertexAttribDivisor(loc, 0);
			glDisableVertexAttribArray(loc);
		}
	}
};


class instance_render_t { // is this a base class of shader_t?

	vector<xform_matrix> inst_xforms;
	int loc;

public:
	instance_render_t(int loc_=-1) : loc(loc_) {}
	void set_loc(int loc_) {loc = loc_;}
	void add_cur_inst();
	void add_inst(xform_matrix const &xf) {inst_xforms.push_back(xf);}
	void draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo=0, int index_type=GL_NONE, void *indices=NULL);
	unsigned size () const {return inst_xforms.size();}
	bool     empty() const {return inst_xforms.empty();}
};


#endif // _SHADERS_H_


