// 3D World
// by Frank Gennari
// Lighting/Lightmap supporting classes
// 1/19/06
#pragma once

#include "3DWorld.h"
#include "trigger.h"

extern int MESH_X_SIZE, MESH_Y_SIZE, MESH_SIZE[3];

#define ADD_LIGHT_CONTRIB(c, C) {C[0] += c[0]; C[1] += c[1]; C[2] += c[2];}

unsigned const FLASHLIGHT_LIGHT_ID = 0;
float const LT_DIR_FALLOFF   = 0.005;
float const LT_DIR_FALLOFF_INV(1.0/LT_DIR_FALLOFF);
float const CTHRESH          = 0.025;
float const SQRT_CTHRESH     = sqrt(CTHRESH);
float const FLASHLIGHT_BW    = 0.02;


class light_grid_base {
protected:
	unsigned get_ix(int x, int y, int z) const {return ((y*MESH_X_SIZE + x)*MESH_SIZE[2] + z);}
	int check_lmap_get_grid_index(point const &p) const;
};


unsigned const lmcell_ltype_off[NUM_LIGHTING_TYPES] = {0, 4, 8, 0}; // sky, global, local, sky cobj accum, dynamic

struct lmcell { // size = 52

	float sc[3], sv, gc[3], gv, lc[3], smoke; // *c[3]: RGB sky, global, local colors
	unsigned char pflow[3]; // flow: x, y, z
	
	lmcell() : sv(0.0), gv(0.0), smoke(0.0) {UNROLL_3X(sc[i_] = gc[i_] = lc[i_] = 0.0; pflow[i_] = 255;)}
	float       *get_offset(int ltype)       {return (sc + lmcell_ltype_off[ltype]);}
	float const *get_offset(int ltype) const {return (sc + lmcell_ltype_off[ltype]);}
	static unsigned get_dsz(int ltype)       {return ((ltype == LIGHTING_LOCAL) ? 3 : 4);}
	void get_final_color(colorRGB &color, float max_indir, float indir_scale=1.0, float extra_ambient=0.0) const;
	void get_final_color_local(colorRGB &color) const;
	void set_outside_colors();
	void mix_lighting_with(lmcell const &lmc, float val);
};


class lmap_manager_t {

	vector<lmcell> vldata_alloc;
	unsigned lm_xsize, lm_ysize, lm_zsize;
	lmcell ***vlmap; // y, x, z (size is determined by {MESH_Y_SIZE, MESH_X_SIZE, MESH_Z_SIZE}

	lmap_manager_t(lmap_manager_t const &) = delete; // forbidden
	void operator=(lmap_manager_t const &) = delete; // forbidden

public:
	bool was_updated;
	cube_t update_bcube;

	lmap_manager_t() : lm_xsize(0), lm_ysize(0), lm_zsize(0), vlmap(NULL), was_updated(0) {update_bcube.set_to_zeros();}
	void clear_cells() {vldata_alloc.clear();} // vlmap matrix headers are not cleared
	bool is_allocated() const {return (vlmap != NULL && !vldata_alloc.empty());}
	size_t size() const {return vldata_alloc.size();}
	bool read_data_from_file(char const *const fn, int ltype);
	bool write_data_to_file(char const *const fn, int ltype) const;
	void clear_lighting_values(int ltype);
	bool is_valid_cell(int x, int y, int z) const;
	lmcell const *get_column(int x, int y) const {return vlmap[y][x];} // Note: no bounds checking
	lmcell *get_column(int x, int y) {return vlmap[y][x];} // Note: no bounds checking
	lmcell &get_lmcell(int x, int y, int z) {return get_column(x, y)[z];} // Note: no bounds checking
	lmcell *get_lmcell_round_down(point const &p);
	lmcell *get_lmcell(point const &p);
	void reset_all(lmcell const &init_lmcell=lmcell());
	template<typename T> void alloc(unsigned nbins, unsigned xsize, unsigned ysize, unsigned zsize, T **nonempty_bins, lmcell const &init_lmcell);
	void init_from(lmap_manager_t const &src);
	void copy_data(lmap_manager_t const &src, float blend_weight=1.0);
};


struct lmcell_local { // size = 12 (must be packed)
	float lc[3];
	lmcell_local() {lc[0] = lc[1] = lc[2] = 0.0;}
	bool is_near_zero(float toler) const {return (lc[0] < toler && lc[1] < toler && lc[2] < toler);}
};

class light_volume_local : public light_grid_base {

	bool changed, compressed;
	unsigned tag_ix;
	float scale; // 0 => disabled
	int bounds[3][2];
	vector<lmcell_local> data;

	unsigned get_num_data() const {return (bounds[0][1] - bounds[0][0])*(bounds[1][1] - bounds[1][0])*(bounds[2][1] - bounds[2][0]);}
	bool read(std::string const &filename);
	bool write(std::string const &filename) const;
	void compress(bool verbose);
public:

	light_volume_local(unsigned tag_ix_) : changed(0), compressed(0), tag_ix(tag_ix_), scale(0.0)
	{bounds[0][0] = bounds[0][1] = bounds[1][0] = bounds[1][1] = bounds[2][0] = bounds[2][1] = 0;}
	void set_bounds(int x1, int x2, int y1, int y2, int z1, int z2);
	void set_scale(float scale_) {changed |= (scale != scale_); scale = scale_;} // changing the scale counts as changed
	bool is_allocated() const {return !data.empty();}
	bool needs_update() const {return (changed     && is_allocated());}
	bool is_active   () const {return (scale > 0.0 && is_allocated());}
	void mark_updated() {changed = 0;}
	void allocate();
	unsigned get_tag_ix() const {return tag_ix;}
	void init(unsigned lvol_ix, float scale_, std::string const &filename);
	void gen_data(unsigned lvol_ix, bool verbose);
	
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
	void write_entry_to_cobj_file(unsigned tag_ix, std::ostream &out) const;
	vector<unsigned> const &get_dlight_ixs_for_tag_ix(unsigned tag_ix) const {
		assert(tag_ix < groups.size());
		return groups[tag_ix].dlight_ixs;
	}
	void create_needed_llvols();
};


struct local_smap_data_t;
class local_smap_manager_t;

class light_source { // size = 116

protected:
	bool dynamic=0, enabled=0, user_placed=0, is_cube_face=0, is_cube_light=0, no_shadows=0;
	unsigned smap_index=0, user_smap_id=0, smap_mgr_id=0, cube_eflags=0, num_dlight_rays=0; // smap_index = index of shadow map texture/data
	float radius=0.0, radius_inv=0.0, r_inner=0.0, bwidth=0.0, near_clip=0.0;
	point pos, pos2; // point/sphere light: use pos; line/cylinder light: use pos and pos2
	vector3d dir;
	colorRGBA color=BLACK;
	cube_t custom_bcube;

	float calc_cylin_end_radius(float falloff=0.0) const;
	local_smap_manager_t &get_smap_mgr() const;

public:
	light_source() {}
	light_source(float sz, point const &p, point const &p2, colorRGBA const &c, bool id=0, vector3d const &d=zero_vector, float bw=1.0, float ri=0.0, bool icf=0, float nc=0.0);
	void mark_is_cube_light(unsigned eflags) {is_cube_light = 1; cube_eflags = eflags;}
	void set_dynamic_state(point const &pos_, vector3d const &dir_, colorRGBA const &color_, bool enabled_) {pos = pos2 = pos_; dir = dir_; color = color_; enabled = enabled_;}
	void set_num_dlight_rays(unsigned num) {num_dlight_rays = num;} // zero = use default
	void add_color(colorRGBA const &c);
	colorRGBA const &get_color() const {return color;}
	float get_radius()           const {return radius;}
	float get_r_inner()          const {return r_inner;} // > 0.0 for sphere light
	float get_near_clip()        const {return near_clip;}
	float get_beamwidth()        const {return bwidth;}
	point const &get_pos()       const {return pos;}
	point const &get_pos2()      const {return pos2;}
	vector3d const &get_dir()    const {return dir;}
	sphere_t get_bsphere()       const {return sphere_t(pos, radius);}
	float get_intensity_at(point const &p, point &updated_lpos) const;
	float get_dir_intensity(vector3d const &obj_dir) const;
	void get_bounds(cube_t &bcube, int bnds[3][2], float sqrt_thresh, bool clip_to_scene_bcube=0, vector3d const &bounds_offset=zero_vector) const;
	void set_custom_bcube(cube_t const c) {custom_bcube = c;}
	cube_t calc_bcube(bool add_pad=0, float sqrt_thresh=0.0, bool clip_to_scene_bcube=0, float falloff=0.0) const;
	cylinder_3dw calc_bounding_cylin( float sqrt_thresh=0.0, bool clip_to_scene_bcube=0, float falloff=0.0) const;
	pos_dir_up calc_pdu(bool dynamic_cobj, bool is_cube_face, float falloff) const;
	unsigned get_cube_eflags() const {return cube_eflags;}
	unsigned get_num_rays()    const {return num_dlight_rays;}
	bool is_visible()     const;
	bool is_directional() const {return (bwidth < 1.0f);}
	bool is_very_directional() const {return ((bwidth + LT_DIR_FALLOFF) < 0.5f);}
	bool is_line_light()  const {return (pos != pos2 && !is_cube_light);} // technically cylinder light
	bool get_is_cube_light() const {return is_cube_light;}
	bool is_dynamic()     const {return dynamic;}
	bool is_neg_light()   const {return (color.R < 0.0 || color.G < 0.0 || color.B < 0.0);}
	bool is_enabled()     const {return enabled;}
	bool is_user_placed() const {return user_placed;}
	bool smap_enabled()   const {return (smap_index != 0 || is_cube_face);}
	bool has_no_shadows() const {return no_shadows;}
	bool is_enabled_spotlight() const {return (is_enabled() && !is_cube_face && !is_cube_light && !is_line_light() && is_very_directional());}
	void set_enabled(bool enabled_) {enabled = enabled_;}
	void disable_shadows() {no_shadows = 1;}
	void shift_by(vector3d const &vd) {pos += vd; pos2 += vd;}
	void pack_to_floatv(float *data) const;
	void combine_with(light_source const &l);
	bool try_merge_into(light_source &ls) const;
	void setup_and_bind_smap_texture(shader_t &s, bool &arr_tex_set) const;
	void write_to_cobj_file(std::ostream &out, bool is_diffuse) const;
	void draw_light_cone(shader_t &shader, float alpha) const;
	bool alloc_shadow_map(bool &matched_smap_id, unsigned sm_size=0);
	void update_shadow_map(bool matched_smap_id, float falloff, bool dynamic_cobj=0, bool outdoor_shadows=0, bool force_update=0);
	bool setup_shadow_map(float falloff, bool dynamic_cobj=0, bool outdoor_shadows=0, bool force_update=0, unsigned sm_size=0);
	void release_smap();
	void invalidate_cached_smap_id(unsigned smap_id) const;
	void assign_smap_id    (unsigned id) {user_smap_id = id;}
	void assign_smap_mgr_id(unsigned id) {smap_mgr_id  = id;}
	bool operator<(light_source const &l) const {return (radius < l.radius);} // compare radius
	bool operator>(light_source const &l) const {return (radius > l.radius);} // compare radius
};


class bind_point_t {

protected:
	bool bound, valid, disabled, dynamic_cobj;
	int bind_cobj;
	point bind_pos;

public:
	bind_point_t() : bound(0), valid(1), disabled(0), dynamic_cobj(0), bind_cobj(-1) {}
	bind_point_t(point const &pos, bool dynamic_=0) : bound(1), valid(1), disabled(0), dynamic_cobj(dynamic_), bind_cobj(-1), bind_pos(pos) {}
	void disable() {disabled = 1;}
	void bind_to_pos(point const &pos, bool dynamic_=0, int bind_cobj_=-1) {bind_pos = pos; bound = 1; dynamic_cobj = dynamic_; bind_cobj = bind_cobj_;}
	bool is_valid();
	point get_updated_bind_pos() const;
	void shift_by(vector3d const &vd) {bind_pos += vd;} // invalidate bind_cobj?
};


class light_source_trig : public light_source, public bind_point_t {

	bool use_smap, outdoor_shadows, dynamic_indir;
	short platform_id;
	unsigned indir_dlight_ix, sm_size; // Note: sm_size of 0 uses default shadow map resolution
	float active_time, inactive_time;
	point last_pos;
	vector3d last_dir;
	multi_trigger_t triggers;
	sensor_t sensor;

	float rot_rate;
	vector3d rot_axis;

public:
	light_source_trig() : use_smap(0), outdoor_shadows(0), dynamic_indir(0), platform_id(-1), indir_dlight_ix(0), sm_size(0),
		active_time(0.0), inactive_time(0.0), last_pos(all_zeros), last_dir(zero_vector), rot_rate(0.0), rot_axis(zero_vector) {}
	light_source_trig(light_source const &ls, bool smap=0, short platform_id_=-1, unsigned lix=0, sensor_t const &cur_sensor=sensor_t(), bool outdoor_shadows_=0, unsigned sm_size_=0)
		: light_source(ls), use_smap(smap), outdoor_shadows(outdoor_shadows_), dynamic_indir(0), platform_id(platform_id_), indir_dlight_ix(lix), sm_size(sm_size_),
		active_time(0.0), inactive_time(0.0), last_pos(pos), last_dir(dir), sensor(cur_sensor), rot_rate(0.0), rot_axis(zero_vector)
	{user_placed = 1; dynamic = (platform_id >= 0); if (is_cube_face) {assert(use_smap);}}
	void add_triggers(multi_trigger_t const &t) {triggers.add_triggers(t);} // deep copy
	void set_rotate(vector3d const &axis, float rotate);
	void enable_dynamic_indir() {dynamic_indir = 1;}
	bool check_activate(point const &p, float radius, int activator);
	void register_activate(bool player_triggered);
	void advance_timestep();
	bool is_enabled() {return (bind_point_t::is_valid() && light_source::is_enabled());}
	void disable() {release_smap(); bind_point_t::disable();}
	bool has_bound_platform() const {return (platform_id >= 0);}
	void shift_by(vector3d const &vd);
	void move_to(point const &new_pos) {shift_by(new_pos - pos);}
	bool is_shadow_map_enabled() const;
	bool check_shadow_map();
	unsigned get_indir_dlight_ix() const {return indir_dlight_ix;}
	bool need_update_indir(); // Note: modifies last_pos/last_dir, not const
	void write_to_cobj_file(std::ostream &out, bool is_diffuse) const;
};


unsigned const MAX_LSRC = 255; // max of 255 lights per bin

class dls_cell {

	unsigned short lsrc[MAX_LSRC+1] = {0};
	unsigned sz;

public:
	dls_cell() : sz(0) {}
	void clear() {sz = 0;}
	
	void add_light(unsigned ix, unsigned char &enabled_flag) {
		if (!enabled_flag) {sz = 0; enabled_flag = 1;} // clear if marked as disabled, then enable
		if (sz < MAX_LSRC) {lsrc[sz++] = ix;}
	}
	void add_light_range(unsigned six, unsigned eix, unsigned char &enabled_flag);
	bool check_add_light(unsigned ix) const;
	size_t size() const {return sz;}
	bool empty()  const {return (sz == 0);}
	unsigned get(unsigned i) const {return lsrc[i];} // no bounds checking
	unsigned short const *get_src_ixs() const {return lsrc;}
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


// from ray_trace.cpp
void check_for_lighting_finished();
void compute_ray_trace_lighting(unsigned ltype, bool verbose);
unsigned add_path_to_lmcs(lmap_manager_t *lmgr, cube_t *bcube, point p1, point const &p2, float weight, colorRGBA const &color, int ltype, bool first_pt);
// from lightmap.cpp
void update_indir_light_tex_range(lmap_manager_t const &lmap, vector<unsigned char> &tex_data,
	unsigned xsize, unsigned y1, unsigned y2, unsigned zsize, float lighting_exponent=1.0, bool local_only=0, bool mt=0);
void indir_light_tex_from_lmap(unsigned &tid, lmap_manager_t const &lmap, vector<unsigned char> &tex_data,
	unsigned xsize, unsigned ysize, unsigned zsize, float lighting_exponent=1.0, bool local_only=0);

