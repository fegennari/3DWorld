// 3D World - 3D Model Rendering Classes
// by Frank Gennari
// 8/17/11

#ifndef _MODEL3D_H_
#define _MODEL3D_H_

#include "3DWorld.h"

using namespace std;

typedef map<string, unsigned> string_map_t;

colorRGB const def_color(0.0, 0.0, 0.0);


struct geom_xform_t {
	vector3d tv;
	float scale;
	bool mirror[3], swap_dim[3][3];

	geom_xform_t() : tv(zero_vector), scale(1.0) {
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(swap_dim[i][i_] = 0;)
			mirror[i] = 0;
		}
	}
	void xform_pos_rm(point &pos) const {
		UNROLL_3X(if (mirror[i_]) pos[i_] = -pos[i_];)
		
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(if (swap_dim[i][i_]) swap(pos[i], pos[i_]);)
		}
	}
	void xform_pos_rms(point &pos) const {
		xform_pos_rm(pos);
		pos *= scale;
	}
	void xform_pos(point &pos) const {
		xform_pos_rms(pos);
		pos += tv;
	}
	void xform_vect(vector<point> &v) const {
		for (vector<point>::iterator i = v.begin(); i != v.end(); ++i) {
			xform_pos(*i);
		}
	}
};


class vntc_vect_t : public vector<vert_norm_tc> {

	unsigned render_vbo, shadow_vbo;

public:
	vntc_vect_t() : render_vbo(0), shadow_vbo(0) {}
	void render(bool is_shadow_pass) const;
	void render_array(bool is_shadow_pass);
	void free_vbos();
};


struct geom_data_t {

	vector<vntc_vect_t> polygons;
	vntc_vect_t triangles;

	void add_poly(vntc_vect_t const &poly);
	void clear() {polygons.clear();}
	void free_context() {triangles.free_vbos();}
	bool empty() const {return polygons.empty();}
	void render_array(bool is_shadow_pass) {triangles.render_array(is_shadow_pass);}
	void render_polygons(bool is_shadow_pass) const;
};


class texture_manager {

protected:
	deque<texture_t> textures;
	string_map_t tex_map; // maps texture filenames to texture indexes

public:
	unsigned create_texture(string const &fn, bool verbose);
	void clear();
	void free_tids();
	void free_textures();
	void ensure_texture_loaded(texture_t &t) const;
	void ensure_tid_loaded(int tid);
	void ensure_tid_bound(int tid);
	void bind_texture(int tid) const;
};


struct material_t {

	colorRGB ka, kd, ks, ke, tf;
	float ns, ni, alpha, tr;
	unsigned illum;
	int a_tid, d_tid, s_tid, alpha_tid, bump_tid;

	// geometry - does this go here or somewhere else?
	geom_data_t geom;

	material_t() : ka(def_color), kd(def_color), ks(def_color), ke(def_color), tf(def_color), ns(1.0), ni(1.0), alpha(1.0), tr(0.0),
		illum(2), a_tid(-1), d_tid(-1), s_tid(-1), alpha_tid(-1), bump_tid(-1) {}
	int get_render_texture() const {return d_tid;}
	bool is_partial_transparent() const {return (alpha < 1.0 || alpha_tid >= 0);}
	void render(texture_manager const &tm, int default_tid, bool is_shadow_pass);
};


class model3d {

	// geometry
	geom_data_t unbound_geom;
	int unbound_tid;
	colorRGBA unbound_color;

	// materials
	deque<material_t> materials;
	string_map_t mat_map; // maps material names to materials indexes
	set<string> undef_materials; // to reduce warning messages

public:
	// textures
	texture_manager &tm;

	model3d(texture_manager &tm_, int def_tid=-1, colorRGBA const &def_c=WHITE)
		: tm(tm_), unbound_tid((def_tid >= 0) ? def_tid : WHITE_TEX), unbound_color(def_c) {}
	unsigned num_materials(void) const {return materials.size();}

	material_t &get_material(int mat_id) {
		assert(mat_id >= 0 && (unsigned)mat_id < materials.size());
		return materials[mat_id];
	}

	// creation and query
	void add_polygon(vntc_vect_t const &poly, int mat_id);
	int get_material_ix(string const &material_name, string const &fn);
	int find_material(string const &material_name);
	void clear();
	void free_context();
	void load_all_used_tids();
	void bind_all_used_tids();
	void render(bool is_shadow_pass); // const?
};


struct model3ds : public deque<model3d> {

	texture_manager tm;

	void clear();
	void free_context();
	void render(bool is_shadow_pass); // const?
};


bool is_poly_convex(vector<point> const &points);

void free_model_context();
void render_models(bool shadow_pass);

bool read_object_file(char *filename, vector<vector<point> > &ppts, geom_xform_t const &xf,
	int def_tid, colorRGBA const &def_c, bool load_model_file, bool verbose);


#endif // _MODEL3D_H_
