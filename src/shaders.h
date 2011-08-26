// 3D World - Shader Class Declaration
// by Frank Gennari
// 8/4/10
#include "3DWorld.h"

using std::string;


unsigned const TEX0_S_ATTR = 0;
unsigned const TEX0_T_ATTR = 1;


class shader_t {

	unsigned program; // active program
	string prepend_string[3]; // vertex=0, fragment=1, geometry=3
	string prog_name_suffix;
	vector<int> attrib_locs;

	string vs_name, fs_name, gs_name;
	int in_prim, out_prim, verts_out;

	int get_attrib_loc(char const *const name) const;
	int attrib_loc_by_ix(unsigned ix) const;

	void filename_split(string const &fname, vector<string> &fns, char sep) const;
	bool load_shader_file(string const &fname, string &data) const;
	unsigned get_shader(string const &name, unsigned type) const;
	void print_shader_info_log(unsigned shader) const;
	void print_program_info_log() const;

public:
	shader_t() : program(0), in_prim(0), out_prim(0), verts_out(0) {}
	~shader_t() {assert(!program); } // end_shader() should have been called

	void set_vert_shader(string const &vs_name_) {vs_name = vs_name_;}
	void set_frag_shader(string const &fs_name_) {fs_name = fs_name_;}
	void set_geom_shader(string const &gs_name_, int in_prim_, int out_prim_, int verts_out_) {
		gs_name = gs_name_; in_prim = in_prim_; out_prim = out_prim_; verts_out = verts_out_;
	}
	bool begin_shader();
	void end_shader();

	int get_uniform_loc(char const *const name) const;
	void set_uniform_float_array(int loc, float const *const val, unsigned num) const;
	void set_uniform_float      (int loc, float val) const;
	void set_uniform_int        (int loc, int val) const;
	void set_uniform_vector3d   (int loc, vector3d const &val) const;
	void set_uniform_color      (int loc, colorRGBA const &val) const;

	void add_uniform_float_array (char const *const name, float const *const val, unsigned num) const;
	void add_uniform_float       (char const *const name, float val) const;
	void add_uniform_int         (char const *const name, int val) const;
	void add_uniform_vector3d    (char const *const name, vector3d const &val) const;
	void add_uniform_color       (char const *const name, colorRGBA const &val) const;
	bool set_uniform_buffer_data (char const *name, float const *data, unsigned size, unsigned &buffer_id) const;

	void register_attrib_name(char const *name, unsigned bind_ix);
	void add_attrib_float_array(unsigned ix, float const *const val, unsigned num) const;
	void add_attrib_float      (unsigned ix, float val) const;
	void add_attrib_int        (unsigned ix, int val) const;

	void setup_enabled_lights(unsigned num=2, unsigned shaders_enabled=3);
	void setup_scene_bounds() const;
	void setup_fog_scale() const;
	void set_prefix(string const &prefix, unsigned shader_type);
	void set_bool_prefix(string const &name, bool val, unsigned shader_type);
	void set_int_prefix(string const &name, int val, unsigned shader_type);
};

