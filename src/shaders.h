// 3D World - Shader Class Declaration
// by Frank Gennari
// 8/4/10
#pragma once

#include "3DWorld.h"

using std::string;

unsigned const TEX0_S_ATTR  = 0;
unsigned const TEX0_T_ATTR  = 1;
unsigned const TANGENT_ATTR = 2;
enum {SHADER_TYPE_VERT=0, SHADER_TYPE_FRAG, SHADER_TYPE_GEOM, SHADER_TYPE_TESC, SHADER_TYPE_TESE, SHADER_TYPE_COMP, NUM_SHADER_TYPES};

#define make_shader_bool_prefix(name, val) ((val) ? ("const bool " name " = true;") : ("const bool " name " = false;"))

bool check_for_tess_shader();


struct gl_light_params_t {

	point pos, eye_space_pos;
	colorRGBA ambient=ALPHA0, diffuse=ALPHA0, specular=ALPHA0;
	float pos_w=-1.0;
	vector3d atten; // {constant, linear, quadratic}

	gl_light_params_t() {set_atten(-1.0, -1.0, -1.0);}
	void set_a (colorRGBA const &a) {ambient = a;};
	void set_ds(colorRGBA const &c) {diffuse = specular = c;}
	void set_atten(float c, float l, float q) {atten.assign(c, l, q);}
	void set_pos(point const &p, float w);
};


struct zi_unsigned_t {unsigned v=0;};

class property_map_t { // for storing user-defined shader properties
	map<string, string> prop_map;
	static string empty_str;
public:
	void add_property   (string const &key, string const &value) {prop_map[key] = value;}
	bool has_property   (string const &key) const {return (prop_map.find(key) != prop_map.end());}
	bool remove_property(string const &key)       {return (prop_map.erase(key) > 0);}
	
	string const &get_property(string const &key) const {
		auto it(prop_map.find(key));
		return ((it == prop_map.end()) ? empty_str : it->second);
	}
	void clear_properties() {prop_map.clear();}
};

// shader user flags
enum {SHADER_FLAG_NO_ALPHA_TEST=0};

class shader_t : public property_map_t {

	unsigned program=0; // active program
	string prepend_string[NUM_SHADER_TYPES]; // vertex=0, fragment=1, geometry=2, tess_control=3, tess_eval=4, compute=5
	string prog_name_prefix;
	vector<int> attrib_locs;
	string shader_names[NUM_SHADER_TYPES];
	colorRGBA last_spec=ALPHA0;
	float last_metalness=0.0;
	
	struct light_loc_t {
		int v[5];
		bool valid;
		light_loc_t() {invalidate();}
		void invalidate() {v[0] = v[1] = v[2] = v[3] = v[4] = -1; valid = 0;}
	};
	light_loc_t light_locs[MAX_SHADER_LIGHTS];
	gl_light_params_t prev_lps[MAX_SHADER_LIGHTS];
	int vnct_locs[4]; // {vertex, normal, color, tex_coord}
	int emission_loc=-1, specular_color_loc=-1, metalness_loc=-1, ref_ix_loc=-1;

	struct subroutine_val_t {
		vector<unsigned> ixs;
		vector<pair<string, unsigned>> name_to_ix;
		void resize(unsigned sz) {assert(ixs.size() == 0 || ixs.size() == sz); ixs.resize(sz, 0);}
		int get_ix_for_name(char const *const name) const;
	};
	typedef map<unsigned, subroutine_val_t> subroutine_map_t;
	subroutine_map_t subroutines;
	unsigned user_flags=0;
	int pm_loc=-1, mvm_loc=-1, mvmi_loc=-1, mvpm_loc=-1, nm_loc=-1; // matrices

	unsigned get_shader(string const &name, unsigned type) const;
	static void print_shader_info_log(unsigned shader);
	void print_program_info_log() const;
	void cache_vnct_locs();
	void cache_matrix_locs();
	void clear_vntc_locs() {vnct_locs[0] = vnct_locs[1] = vnct_locs[2] = vnct_locs[3] = -1;}

public:
	shader_t() {clear_vntc_locs();}
	//~shader_t() {assert(!program);} // end_shader() should have been called (but not for cached global variables)
	unsigned get_program() const {return program;} // semi-private, for internal use as map key in vao_cache_t
	void get_program_binary(vector<unsigned char> &binary_data, GLenum &binary_format) const;
	void set_program_binary(vector<unsigned char> const &binary_data, GLenum const binary_format);

	void set_vert_shader(string const &vs_name_) {shader_names[0] = vs_name_;}
	void set_frag_shader(string const &fs_name_) {shader_names[1] = fs_name_;}
	void set_geom_shader(string const &gs_name_) {shader_names[2] = gs_name_;}
	void set_tess_control_shader(string const &tcs_name_) {shader_names[3] = tcs_name_;}
	void set_tess_eval_shader   (string const &tes_name_) {shader_names[4] = tes_name_;}
	void set_comp_shader(string const &cs_name_) {shader_names[5] = cs_name_;}
	bool has_tess_shader()   const {return !shader_names[4].empty();}
	bool is_compute_shader() const {return !shader_names[5].empty();}

	bool is_setup() const {return (program > 0);}
	void enable();
	void disable();
	void clear();
	void make_current();
	bool begin_shader(bool do_enable=1);
	void end_shader();
	void begin_color_only_shader();
	void begin_color_only_shader(colorRGBA const &color);
	void begin_shadow_map_shader(bool use_alpha_mask=0, bool enable_xlate_scale=0);
	void begin_simple_textured_shader(float min_alpha=0.0, bool include_2_lights=0, bool use_texgen=0, colorRGBA const *const color=NULL);
	void begin_untextured_lit_glcolor_shader();

	void enable_vnct_atribs(bool va, bool tca, bool na, bool ca) const;
	void set_vertex_ptr(unsigned stride, void const *const ptr) const;
	void set_normal_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_color4_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_tcoord_ptr(unsigned stride, void const *const ptr, bool compressed) const;
	void set_cur_color(colorRGBA const &color) const;
	void set_cur_normal(vector3d const &normal) const;

	int get_uniform_loc(char const *const name) const;
	void ensure_uniform_loc(int &loc, char const *const name) const {if (loc < 0) {loc = get_uniform_loc(name);}}
	static bool set_uniform_float_array(int loc, float const *const val, unsigned num);
	static bool set_uniform_float      (int loc, float val);
	static bool set_uniform_int        (int loc, int val);
	static bool set_uniform_uint       (int loc, unsigned val);
	static bool set_uniform_handle     (int loc, GLuint64 val);
	static bool set_uniform_vector2d   (int loc, vector2d const &val);
	static bool set_uniform_vector3d   (int loc, vector3d const &val);
	static bool set_uniform_vector4d   (int loc, vector4d const &val);
	static bool set_uniform_color      (int loc, colorRGBA const &val);
	static bool set_uniform_color      (int loc, colorRGB  const &val);
	static bool set_uniform_matrix_3x3 (int loc, float const *const m, bool transpose, unsigned num=1);
	static bool set_uniform_matrix_4x4 (int loc, float const *const m, bool transpose, unsigned num=1);

	bool add_uniform_float_array (char const *const name, float const *const val, unsigned num) const;
	bool add_uniform_float       (char const *const name, float val) const;
	bool add_uniform_int         (char const *const name, int val) const;
	bool add_uniform_uint        (char const *const name, unsigned val) const;
	bool add_uniform_handle      (char const *const name, GLuint64 val) const;
	bool add_uniform_vector2d    (char const *const name, vector2d  const &val) const;
	bool add_uniform_vector3d    (char const *const name, vector3d  const &val) const;
	bool add_uniform_vector4d    (char const *const name, vector4d  const &val) const;
	bool add_uniform_color       (char const *const name, colorRGBA const &val) const;
	bool add_uniform_color       (char const *const name, colorRGB  const &val) const;
	bool add_uniform_matrix_4x4  (char const *const name, float const *m, bool transpose, unsigned num=1) const;

	unsigned get_subroutine_index(int shader_type, char const *const name) const;
	unsigned get_subroutine_uniform_loc(int shader_type, char const *const name) const;
	void set_all_subroutines(int shader_type, unsigned count, char const *const *const uniforms, char const *const *const bindings);
	void set_subroutines(int shader_type, unsigned count, unsigned const *const indices);
	void set_subroutines(int shader_type, vector<unsigned> const &indices) {set_subroutines(shader_type, indices.size(), &indices.front());}
	void set_subroutine (int shader_type, unsigned index); // single subroutine
	void set_subroutine (int shader_type, char const *const name) {set_subroutine(shader_type, get_subroutine_index(shader_type, name));}
	void reset_subroutine(int shader_type, char const *const uniform, char const *const binding); // one of multiple
	void restore_subroutines();

	int attrib_loc_by_ix(unsigned ix, bool allow_fail=0) const;
	int get_attrib_loc(char const *const name, bool allow_fail=0) const;
	void register_attrib_name(char const *const name, unsigned bind_ix);
	bool set_attrib_float_array(int loc, float const *const val, unsigned num) const;
	bool set_attrib_int_array  (int loc, int   const *const val, unsigned num) const;
	bool set_attrib_float      (int loc, float val) const;
	bool set_attrib_int        (int loc, int val) const;
	bool add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const;
	bool add_attrib_int_array  (unsigned ix, int   const *const val, unsigned num) const;
	bool add_attrib_float      (unsigned ix, float val) const;
	bool add_attrib_int        (unsigned ix, int val) const;

	void setup_enabled_lights(unsigned num=2, unsigned shaders_enabled=3);
	void upload_light_source(unsigned light_id, unsigned field_filt=0xFF);
	void upload_light_sources_range(unsigned start, unsigned end);
	void upload_all_light_sources() {upload_light_sources_range(0, MAX_SHADER_LIGHTS);}

	void upload_pjm();
	void upload_mvm();

	void setup_scene_bounds() const;
	void setup_scene_bounds_from_bcube(cube_t const &bcube) const;
	void setup_fog_scale() const;
	void check_for_fog_disabled();

	void set_prefix_str(string const &prefix, unsigned shader_type) {set_prefix(prefix.c_str(), shader_type);}
	void set_prefix(char const *const prefix, unsigned shader_type);
	void set_prefixes_str(string const &prefix, unsigned shaders_enabled=3) {set_prefixes(prefix.c_str(), shaders_enabled);}
	void set_prefixes(char const *const prefix, unsigned shaders_enabled=3);
	void set_int_prefix(char const *const name, int val, unsigned shader_type);

	void set_color_e(colorRGBA const &color);
	void clear_color_e() {set_color_e(BLACK);}
	void set_black_diffuse_emissive_color(colorRGBA const &color);
	void set_specular_color(colorRGB const &specular, float shininess);
	void set_specular(float spec, float shine) {set_specular_color(colorRGB(spec, spec, spec), shine);}
	void set_specular(float spec, float shine, float metalness) {set_specular(spec, shine); set_metalness(metalness);}
	void set_metalness(float metalness);
	void set_refract_ix(float ref_ix);
	void set_material(base_mat_t const &mat);
	void clear_specular() {set_specular(0.0, 0.0);}
	void clear_specular_and_metalness() {clear_specular(); set_metalness(0.0);}

	unsigned get_user_flags() const {return user_flags;}
	bool get_user_flag(unsigned fbit) const {assert(fbit < 32); return (user_flags & (1<<fbit));}
	void set_user_flags(unsigned f) {user_flags = f;}
	void set_user_flag(unsigned fbit, bool val=1) {assert(fbit < 32); if (val) {user_flags |= (1<<fbit);} else {user_flags &= ~(1<<fbit);}}
};


// special shader used with volume particle clouds (nebulas, explosions, teleporters) that caches uniform locs
class vpc_shader_t : public shader_t { // move somewhere else?
public:
	int ns_loc=-1, c1i_loc=-1, c1o_loc=-1, c2i_loc=-1, c2o_loc=-1, c3i_loc=-1, c3o_loc=-1, rad_loc=-1, rs_loc=-1, off_loc=-1, vd_loc=-1, as_loc=-1, nexp_loc=-1;
	void cache_locs();
	void set_all_colors(colorRGBA const &c1i, colorRGBA const &c1o, colorRGBA const &c2i, colorRGBA const &c2o, colorRGBA const &c3i, colorRGBA const &c3o);
};


class compute_shader_base_t : public shader_t {
protected:
	unsigned xsize, ysize, xsize_req, ysize_req; // actual and requested sizes
	bool is_running=0;

	bool setup_target_texture(unsigned &tid, bool is_R32F) const;
public:
	compute_shader_base_t(unsigned xsize_, unsigned ysize_) : xsize(xsize_), ysize(ysize_), xsize_req(xsize), ysize_req(ysize) {assert(xsize > 0 && ysize > 0);}
	bool get_is_running() const {return is_running;}
	unsigned get_xsize() const {return xsize;}
	unsigned get_ysize() const {return ysize;}
	unsigned get_xsize_req() const {return xsize_req;}
	unsigned get_ysize_req() const {return ysize_req;}
};

// "fake" compute shader implemented as a fragment shader
class compute_shader_t : public compute_shader_base_t {

	unsigned fbo_id=0, pbo=0;
	string frag_shader_str;

	unsigned get_pbo_size() const {return xsize*ysize*sizeof(float);}
	void draw_geom() const;
	void unset_fbo(bool keep_fbo_for_reuse);
	void read_pixels(vector<float> &vals, bool is_last=1);
public:
	compute_shader_t(string const &fstr, unsigned xsize_, unsigned ysize_) : compute_shader_base_t(xsize_, ysize_), frag_shader_str(fstr) {}
	void begin();
	void end_shader();
	void setup_and_run(unsigned &tid, bool is_R32F, bool is_first=1, bool is_last=1);
	void run(unsigned &tid);
	void prep_for_read_pixels(bool is_first=1);
	void read_float_vals(vector<float> &vals, bool is_last=1, bool keep_fbo_for_reuse=0);
	void gen_matrix_RGBA8(vector<float> &vals, unsigned &tid, bool is_first=1, bool is_last=1, bool keep_fbo_for_reuse=0);
	void gen_matrix_R32F(vector<float> &vals, unsigned &tid, bool is_first=1, bool is_last=1, bool keep_fbo_for_reuse=0);
	void set_comp_prefix(char const *const prefix) {set_prefix(prefix, 1);} // FS
};

// "real" compute shader
class compute_shader_comp_t : public compute_shader_base_t {

	unsigned zsize, zsize_req, block_sz_x, block_sz_y, block_sz_z;
	string comp_shader_str;
public:
	compute_shader_comp_t(string const &cstr, unsigned xsize_, unsigned ysize_, unsigned zsize_=1, unsigned bsx=2, unsigned bsy=2, unsigned bsz=1);
	bool is_3d() const {return (zsize > 1);}
	void begin();
	void setup_and_run(unsigned &tid, bool is_R32F, bool is_first=1, bool is_last=1);
	void prep_for_read_pixels(bool is_first=1) {} // does nothing
	void read_float_vals(vector<float> &vals, bool is_last=1, bool keep_fbo_for_reuse=0);
	void gen_matrix_R32F(vector<float> &vals, unsigned &tid, bool is_first=1, bool is_last=1);
	void set_comp_prefix(char const *const prefix) {set_prefix(prefix, 5);} // CS
};


template<unsigned M, unsigned N> struct shader_float_matrix_uploader {
	static void enable(int start_loc, int divisor, float const *const data=NULL);
	static void disable(int start_loc);
};


class text_drawer {
	shader_t s;
	vector<vert_tc_t> verts;
	colorRGBA cur_color=ALPHA0;
public:
	void begin_draw(colorRGBA const *const color=nullptr);
	void end_draw();
	static void bind_font_texture();
	void set_color(colorRGBA const &color);
	void flush();
	void add_text(string const &text, point const &pos, float tsize, vector3d const &column_dir=plus_x, vector3d const &row_dir=plus_y, colorRGBA const *const color=nullptr);
	void add_text(colorRGBA const &color, float x, float y, float z, char const *text, float tsize) {
		add_text(text, point(x, y, z), 0.8*tsize, plus_x, plus_y, &color);
	}
};

// add some special texture ID enums, used in buildings
unsigned const FONT_TEXTURE_ID = (1<<16); // some large number that will never be a valid texture ID
unsigned const REFLECTION_TEXTURE_ID = FONT_TEXTURE_ID + 1;
unsigned const ABST_ART_TEXTURE_ID   = FONT_TEXTURE_ID + 2;

struct tile_blend_tex_data_t {
	unsigned tid_tinput=0, tid_lut=0, context_count=0;
	vector3d colorSpaceVector1, colorSpaceVector2, colorSpaceVector3, colorSpaceOrigin;

	bool textures_valid() const {return (tid_tinput > 0 && tid_lut > 0);}
	void create_textures(texture_t const &texture);
	void ensure_textures(unsigned tid);
	void bind_shader(shader_t &s) const;
	void clear_context();
};

void set_one_texture(shader_t &s, unsigned tid, unsigned tu_id, const char *const name);
void setup_shader_underwater_atten(shader_t &s, float atten_scale, float mud_amt=0.0, float algae_amt=0.0);

