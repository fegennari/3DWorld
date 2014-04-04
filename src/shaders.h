// 3D World - Shader Class Declaration
// by Frank Gennari
// 8/4/10

#ifndef _SHADERS_H_
#define _SHADERS_H_

#include "3DWorld.h"

using std::string;


unsigned const TEX0_S_ATTR       = 0;
unsigned const TEX0_T_ATTR       = 1;
unsigned const NUM_SHADER_TYPES  = 5;


struct gl_light_params_t {

	point pos, eye_space_pos;
	colorRGBA ambient, diffuse, specular;
	float pos_w;
	float const_atten;
	float linear_atten;
	float quad_atten;

	gl_light_params_t() : pos(all_zeros), eye_space_pos(pos), ambient(ALPHA0), diffuse(ALPHA0), specular(ALPHA0), pos_w(-1.0) {set_atten(-1.0, -1.0, -1.0);}
	void set_a (colorRGBA const &a) {ambient = a;};
	void set_ds(colorRGBA const &c) {diffuse = specular = c;}
	void set_atten(float c, float l, float q) {const_atten = c; linear_atten = l; quad_atten = q;}
	void set_pos(point const &p, float w);
};


class shader_t {

	unsigned program; // active program
	string prepend_string[NUM_SHADER_TYPES]; // vertex=0, fragment=1, geometry=2, tess_control=3, tess_eval=4
	string prog_name_prefix;
	vector<int> attrib_locs;
	string shader_names[NUM_SHADER_TYPES];
	colorRGBA last_spec;

	struct light_loc_t {
		int v[7];
		bool valid;
		light_loc_t() : valid(0) {}
	};
	light_loc_t light_locs[MAX_SHADER_LIGHTS];
	gl_light_params_t prev_lps[MAX_SHADER_LIGHTS];
	int vnct_locs[4]; // {vertex, normal, color, tex_coord}

	unsigned get_shader(string const &name, unsigned type) const;
	static void print_shader_info_log(unsigned shader);
	void print_program_info_log() const;
	void cache_vnct_locs();

public:
	shader_t() : program(0), last_spec(ALPHA0) {}
	//~shader_t() {assert(!program);} // end_shader() should have been called (but not for cached global variables)

	void set_vert_shader(string const &vs_name_) {shader_names[0] = vs_name_;}
	void set_frag_shader(string const &fs_name_) {shader_names[1] = fs_name_;}
	void set_geom_shader(string const &gs_name_) {shader_names[2] = gs_name_;}
	void set_tess_control_shader(string const &tcs_name_) {shader_names[3] = tcs_name_;}
	void set_tess_eval_shader   (string const &tes_name_) {shader_names[4] = tes_name_;}

	bool is_setup() const {return (program > 0);}
	void enable();
	void disable();
	bool begin_shader(bool do_enable=1);
	void end_shader();
	void begin_color_only_shader();
	void begin_color_only_shader(colorRGBA const &color);
	void begin_simple_textured_shader(float min_alpha=0.0, bool include_2_lights=0, bool use_texgen=0, colorRGBA const *const color=NULL);
	void begin_untextured_lit_glcolor_shader();

	void enable_vnct_atribs(bool va, bool tca, bool na, bool ca) const;
	void set_vertex_ptr(unsigned stride, void const *const ptr) const;
	void set_normal_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_color4_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_tcoord_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_cur_color(colorRGBA const &color) const;

	int get_uniform_loc(char const *const name) const;
	static bool set_uniform_float_array(int loc, float const *const val, unsigned num);
	static bool set_uniform_float      (int loc, float val);
	static bool set_uniform_int        (int loc, int val);
	static bool set_uniform_vector2d   (int loc, vector2d const &val);
	static bool set_uniform_vector3d   (int loc, vector3d const &val);
	static bool set_uniform_vector4d   (int loc, vector4d const &val);
	static bool set_uniform_color      (int loc, colorRGBA const &val);
	static bool set_uniform_color      (int loc, colorRGB  const &val);
	static bool set_uniform_matrid_4x4 (int loc, float *m, bool transpose);

	bool add_uniform_float_array (char const *const name, float const *const val, unsigned num) const;
	bool add_uniform_float       (char const *const name, float val) const;
	bool add_uniform_int         (char const *const name, int val) const;
	bool add_uniform_vector2d    (char const *const name, vector2d const &val) const;
	bool add_uniform_vector3d    (char const *const name, vector3d const &val) const;
	bool add_uniform_color       (char const *const name, colorRGBA const &val) const;
	bool add_uniform_color       (char const *const name, colorRGB  const &val) const;
	bool add_uniform_matrix_4x4  (char const *const name, float *m, bool transpose) const;
	bool set_uniform_buffer_data (char const *name, float const *data, unsigned size, unsigned &buffer_id) const;

	int attrib_loc_by_ix(unsigned ix) const;
	int get_attrib_loc(char const *const name, bool allow_fail=0) const;
	void register_attrib_name(char const *const name, unsigned bind_ix);
	bool set_attrib_float_array(int loc, float const *const val, unsigned num) const;
	bool set_attrib_float      (int loc, float val) const;
	bool set_attrib_int        (int loc, int val) const;
	bool add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const;
	bool add_attrib_float      (unsigned ix, float val) const;
	bool add_attrib_int        (unsigned ix, int val) const;

	void setup_enabled_lights(unsigned num=2, unsigned shaders_enabled=3);
	void upload_light_source(unsigned light_id, unsigned field_filt=0xFF);
	void upload_light_sources_range(unsigned start, unsigned end);
	void upload_all_light_sources() {upload_light_sources_range(0, MAX_SHADER_LIGHTS);}
	void setup_scene_bounds() const;
	void setup_fog_scale() const;
	void check_for_fog_disabled();
	void set_prefix_str(string const &prefix, unsigned shader_type) {set_prefix(prefix.c_str(), shader_type);}
	void set_prefix(char const *const prefix, unsigned shader_type);
	void set_bool_prefix(char const *const name, bool val, unsigned shader_type);
	void set_bool_prefixes(char const *const name, bool val, unsigned shaders_enabled=3);
	void set_int_prefix(char const *const name, int val, unsigned shader_type);

	void set_color_e(colorRGBA const &color);
	void clear_color_e() {set_color_e(BLACK);}
	void set_specular_color(colorRGB const &specular, float shininess);
	void set_specular(float spec, float shine) {set_specular_color(colorRGB(spec, spec, spec), shine);}
};


template<unsigned M, unsigned N> struct shader_float_matrix_uploader {

	static void enable(int start_loc, int divisor, float const *const data=NULL) {
		assert(start_loc >= 0 && divisor >= 0);

		for (unsigned i = 0; i < N; ++i) {
			int const loc(start_loc + i);
			glEnableVertexAttribArray(loc);
			glVertexAttribPointer(loc, M, GL_FLOAT, GL_FALSE, M*N*sizeof(float), (const void *)(data + M*i));
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


#endif // _SHADERS_H_


