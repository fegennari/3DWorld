// 3D World
// by Frank Gennari
// Lighting/Lightmap supporting classes
// 1/19/06
#ifndef _LIGHTMAP_H_
#define _LIGHTMAP_H_

#include "3DWorld.h"
#include "trigger.h"

extern int MESH_X_SIZE, MESH_Y_SIZE, MESH_SIZE[3];

#define ADD_LIGHT_CONTRIB(c, C) {C[0] += c[0]; C[1] += c[1]; C[2] += c[2];}

unsigned const FLASHLIGHT_LIGHT_ID = 0;

float const LT_DIR_FALLOFF   = 0.005;
float const LT_DIR_FALLOFF_INV(1.0/LT_DIR_FALLOFF);
float const CTHRESH          = 0.025;
float const SQRT_CTHRESH     = sqrt(CTHRESH);


struct normal_cell { // size = 24, unused

	vector3d n[2]; // {negative, positive}
	normal_cell() {UNROLL_3X(n[0][i_] = n[1][i_] = 0.0;)}

	void add_normal(vector3d const &N, float weight=1.0) { // normal should be normalized
		UNROLL_3X(n[N[i_] > 0.0][i_] += weight*N[i_];)
	}
	void normalize() {
		float const mag(sqrt(n[0].mag_sq() + n[1].mag_sq()));
		n[0] /= mag; n[1] /= mag;
	}
	float dot_product(vector3d const &v) const {
		float dp(0.0);
		UNROLL_3X(dp += n[v[i_] > 0.0][i_]*v[i_];)
		assert(dp >= 0.0);
		return dp;
	}
};


struct lmcell { // size = 52

	float sc[3], sv, gc[3], gv, lc[3], smoke; // *c[3]: RGB sky, global, local colors
	unsigned char pflow[3]; // flow: x, y, z
	
	lmcell() : sv(0.0), gv(0.0), smoke(0.0) {UNROLL_3X(sc[i_] = gc[i_] = lc[i_] = 0.0; pflow[i_] = 255;)}
	float       *get_offset(int ltype)       {return (sc + 4*ltype);}
	float const *get_offset(int ltype) const {return (sc + 4*ltype);}
	static unsigned get_dsz(int ltype)       {return ((ltype == LIGHTING_LOCAL) ? 3 : 4);}
	void get_final_color(colorRGB &color, float max_indir, float indir_scale=1.0, float extra_ambient=0.0) const;
	void set_outside_colors();
	void mix_lighting_with(lmcell const &lmc, float val);
};


class lmap_manager_t {

	vector<lmcell> vldata_alloc;
	unsigned lm_zsize;
	lmcell ***vlmap; // y, x, z (size is determined by {MESH_Y_SIZE, MESH_X_SIZE, MESH_Z_SIZE}

	lmap_manager_t(lmap_manager_t const &); // forbidden
	void operator=(lmap_manager_t const &); // forbidden

public:
	bool was_updated;

	lmap_manager_t() : lm_zsize(0), vlmap(NULL), was_updated(0) {}
	void clear_cells() {vldata_alloc.clear();} // vlmap matrix headers are not cleared
	bool is_allocated() const {return (vlmap != NULL && !vldata_alloc.empty());}
	size_t size() const {return vldata_alloc.size();}
	bool read_data_from_file(char const *const fn, int ltype);
	bool write_data_to_file(char const *const fn, int ltype) const;
	void clear_lighting_values(int ltype);
	bool is_valid_cell(int x, int y, int z) const;
	lmcell *get_column(int x, int y) {return vlmap[y][x];} // Note: no bounds checking
	lmcell &get_lmcell(int x, int y, int z) {return get_column(x, y)[z];} // Note: no bounds checking
	lmcell *get_lmcell(point const &p);
	template<typename T> void alloc(unsigned nbins, unsigned zsize, T **nonempty_bins, lmcell const &init_lmcell);
	void init_from(lmap_manager_t const &src);
	void copy_data(lmap_manager_t const &src, float blend_weight=1.0);
};


struct lmcell_local { // size = 12 (must be packed)
	float lc[3];
	lmcell_local() {lc[0] = lc[1] = lc[2] = 0.0;}
	bool is_near_zero(float toler) const {return (lc[0] < toler && lc[1] < toler && lc[2] < toler);}
};

class light_volume_local {

	bool changed, compressed;
	unsigned tag_ix;
	float scale; // 0 => disabled
	int bounds[3][2];
	vector<lmcell_local> data;

	unsigned get_num_data() const {return (bounds[0][1] - bounds[0][0])*(bounds[1][1] - bounds[1][0])*(bounds[2][1] - bounds[2][0]);}
	bool read(std::string const &filename);
	bool write(std::string const &filename) const;
	void compress();
	unsigned get_ix(int x, int y, int z) const {return ((y*MESH_X_SIZE + x)*MESH_SIZE[2] + z);}
public:

	light_volume_local(unsigned tag_ix_) : changed(0), compressed(0), tag_ix(tag_ix_), scale(0.0) {}
	void set_bounds(int x1, int x2, int y1, int y2, int z1, int z2);
	void set_scale(float scale_) {changed |= (scale != scale_); scale = scale_;} // changing the scale counts as changed
	bool is_allocated() const {return !data.empty();}
	bool needs_update() const {return (changed     && is_allocated());}
	bool is_active   () const {return (scale > 0.0 && is_allocated());}
	void mark_updated() {changed = 0;}
	void allocate();
	unsigned get_tag_ix() const {return tag_ix;}
	void init(unsigned lvol_ix, float scale_, std::string const &filename);
	
	void reset_to_zero() {
		if (!is_allocated()) return;
		data.clear(); allocate(); // clear + resize should re-construct the cells to all zeros
		changed = 1;
	}
	bool check_xy_bounds(int x, int y) const {return (x >= bounds[0][0] && x < bounds[0][1] && y >= bounds[1][0] && y < bounds[1][1]);}
	void add_lighting(colorRGB &color, int x, int y, int z) const;
	void add_color(point const &p, colorRGBA const &color);
};

typedef vector<std::unique_ptr<light_volume_local>> llv_vect;


class tag_ix_map {

	map<std::string, unsigned> name_to_ix;
	unsigned next_ix;

public:
	tag_ix_map() : next_ix(1) {} // ix starts at 1, 0 is a special value for "none"
	unsigned get_ix_for_name(std::string const &name);
};

class indir_dlight_group_manager_t : public tag_ix_map {

	struct group_t {
		int llvol_ix;
		float scale;
		std::string filename;
		vector<unsigned> dlight_ixs; // Note: dynamic lights should all share the same trigger
		group_t(float scale_=1.0) : llvol_ix(-1), scale(scale_) {}
	};
	vector<group_t> groups;
public:
	unsigned get_ix_for_name(std::string const &name, float scale=1.0);
	void add_dlight_ix_for_tag_ix(unsigned tag_ix, unsigned dlight_ix);
	vector<unsigned> const &get_dlight_ixs_for_tag_ix(unsigned tag_ix) const {
		assert(tag_ix < groups.size());
		return groups[tag_ix].dlight_ixs;
	}
	void create_needed_llvols();
};


struct local_smap_data_t;

class light_source { // size = 68

protected:
	bool dynamic, enabled, user_placed;
	unsigned smap_index; // index of shadow map texture/data
	float radius, radius_inv, r_inner, bwidth;
	point pos, pos2; // point/sphere light: use pos; line/cylinder light: use pos and pos2
	vector3d dir;
	colorRGBA color;

	float calc_cylin_end_radius() const;

public:
	light_source() : enabled(0), user_placed(0), smap_index(0) {}
	light_source(float sz, point const &p, point const &p2, colorRGBA const &c, bool id, vector3d const &d=zero_vector, float bw=1.0, float ri=0.0);
	void set_dynamic_state(point const &pos_, vector3d const &dir_, bool enabled_) {pos = pos2 = pos_; dir = dir_; enabled = enabled_;}
	void add_color(colorRGBA const &c);
	colorRGBA const &get_color() const {return color;}
	float get_radius()           const {return radius;}
	float get_r_inner()          const {return r_inner;} // > 0.0 for sphere light
	point const &get_pos()       const {return pos;}
	point const &get_pos2()      const {return pos2;}
	float get_intensity_at(point const &p, point &updated_lpos) const;
	float get_dir_intensity(vector3d const &obj_dir) const;
	void get_bounds(cube_t &bcube, int bnds[3][2], float sqrt_thresh, vector3d const &bounds_offset=zero_vector) const;
	cube_t calc_bcube(float sqrt_thresh=0.0) const;
	cylinder_3dw calc_bounding_cylin(float sqrt_thresh=0.0) const;
	pos_dir_up calc_pdu() const;
	bool is_visible()     const;
	bool is_directional() const {return (bwidth < 1.0);}
	bool is_very_directional() const {return ((bwidth + LT_DIR_FALLOFF) < 0.5);}
	bool is_line_light()  const {return (pos != pos2);} // technically cylinder light
	bool is_dynamic()     const {return dynamic;}
	bool is_neg_light()   const {return (color.R < 0.0 || color.G < 0.0 || color.B < 0.0);}
	bool is_enabled()     const {return enabled;}
	bool is_user_placed() const {return user_placed;}
	bool smap_enabled()   const {return (smap_index != 0);}
	void set_enabled(bool enabled_) {enabled = enabled_;}
	void shift_by(vector3d const &vd) {pos += vd; pos2 += vd;}
	void pack_to_floatv(float *data) const;
	void combine_with(light_source const &l);
	bool try_merge_into(light_source &ls) const;
	void setup_and_bind_smap_texture(shader_t &s, bool &arr_tex_set) const;
	bool operator<(light_source const &l) const {return (radius < l.radius);} // compare radius
	bool operator>(light_source const &l) const {return (radius > l.radius);} // compare radius
};


class bind_point_t {

protected:
	bool bound, valid;
	int bind_cobj;
	point bind_pos;

public:
	bind_point_t() : bound(0), valid(1), bind_cobj(-1) {}
	bind_point_t(point const &pos) : bound(1), valid(1), bind_cobj(-1), bind_pos(pos) {}
	void bind_to_pos(point const &pos) {bind_pos = pos; bound = 1;}
	bool is_valid();
	void shift_by(vector3d const &vd) {bind_pos += vd;} // invalidate bind_cobj?
};


class light_source_trig : public light_source, public bind_point_t {

	bool use_smap;
	short platform_id;
	unsigned indir_dlight_ix;
	float active_time, inactive_time;
	multi_trigger_t triggers;

public:
	light_source_trig() : use_smap(0), platform_id(-1), indir_dlight_ix(0), active_time(0.0), inactive_time(0.0) {}
	light_source_trig(light_source const &ls, bool smap=0, short platform_id_=-1, unsigned lix=0)
		: light_source(ls), use_smap(smap), platform_id(platform_id_), indir_dlight_ix(lix), active_time(0.0), inactive_time(0.0) {user_placed = 1; dynamic = (platform_id >= 0);}
	void add_triggers(multi_trigger_t const &t) {triggers.add_triggers(t);} // deep copy
	bool check_activate(point const &p, float radius, int activator);
	void advance_timestep();
	bool is_enabled() {return (bind_point_t::is_valid() && light_source::is_enabled());}
	void shift_by(vector3d const &vd);
	bool is_shadow_map_enabled() const;
	bool check_shadow_map();
	void release_smap();
	unsigned get_indir_dlight_ix() const {return indir_dlight_ix;}
};


class dls_cell {

	float z1, z2;
	vector<unsigned short> lsrc;

public:
	dls_cell() : z1(FAR_DISTANCE), z2(-FAR_DISTANCE) {}

	void clear();
	void add_light(unsigned ix, float zmin, float zmax);
	bool check_add_light(unsigned ix) const;
	size_t size() const {return lsrc.size();}
	bool empty()  const {return lsrc.empty();}
	unsigned get(unsigned i) const {return lsrc[i];} // no bounds checking
	bool check_z(float z)    const {return (!empty() && z >= z1 && z <= z2);}
	bool check_z_range(float zlo, float zhi) const {return (!empty() && zlo <= z2 && zhi >= z1);}
	vector<unsigned short> const &get_src_ixs() const {return lsrc;}
};


struct cube_light_src {

	cube_t bounds;
	colorRGB color;
	float intensity;
	unsigned num_rays, disabled_edges;

	cube_light_src() : color(BLACK), intensity(0.0), num_rays(0), disabled_edges(0) {}
};


class cube_light_src_vect : public vector<cube_light_src> {
public:
	bool ray_intersects_any(point const &start_pt, point const &end_pt) const;
};


void check_for_lighting_finished();
void compute_ray_trace_lighting(unsigned ltype);


#endif

