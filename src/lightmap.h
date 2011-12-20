// 3D World
// by Frank Gennari
// Lighting/Lightmap supporting classes
// 1/19/06
#ifndef _LIGHTMAP_H_
#define _LIGHTMAP_H_

#include "3DWorld.h"

typedef short CELL_LOC_T;

point const init_point(FAR_CLIP, FAR_CLIP, FAR_CLIP);
extern int MESH_SIZE[3];

#define ADD_LIGHT_CONTRIB(c, C) {C[0] += c[0]; C[1] += c[1]; C[2] += c[2];}


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
	void pack(unsigned char *data, unsigned &pos) const {
		assert(data);
		for (unsigned d = 0; d < 2; ++d) {
			UNROLL_3X(data[pos++] = (unsigned char)(255.0*CLIP_TO_01((d ? 1.0f : -1.0f)*n[d][i_]));)
		}
	}
	void unpack(unsigned char const *data, unsigned &pos) {
		assert(data);
		for (unsigned d = 0; d < 2; ++d) {
			UNROLL_3X(n[d][i_] = (d ? 1.0 : -1.0)*data[pos++]/255.0;)
		}
	}
};


struct lmcell { // size = 56

	float sc[3], sv, gc[3], gv, lc[3], smoke; // *c[3]: RGB sky, global, local colors
	unsigned char lflow[3], pflow[3]; // flow: x, y, z
	//normal_cell n;
	
	lmcell() : sv(0.0), gv(0.0), smoke(0.0) {
		UNROLL_3X(sc[i_] = gc[i_] = lc[i_] = 0.0; lflow[i_] = pflow[i_] = 255;)
	}
	float       *get_offset(int ltype)       {return (sc + 4*ltype);}
	float const *get_offset(int ltype) const {return (sc + 4*ltype);}
	static unsigned get_dsz(int ltype)       {return ((ltype == LIGHTING_LOCAL) ? 3 : 4);}
	void get_final_color(colorRGB &color, float max_indir) const;
	void set_outside_colors();
};


class lmap_manager_t {

	vector<lmcell> vldata_alloc;
	unsigned lm_zsize;

public:
	lmcell ***vlmap; // y, x, z

	lmap_manager_t() : lm_zsize(0), vlmap(NULL) {}
	void clear() {vldata_alloc.clear();} // reset vlmap to NULL?
	bool is_allocated() const {return (vlmap != NULL);}
	size_t size() const {return vldata_alloc.size();}
	bool read_data_from_file(char const *const fn, int ltype);
	bool write_data_to_file(char const *const fn, int ltype) const;
	void clear_lighting_values(int ltype);

	inline bool is_valid_cell(int x, int y, int z) const {
		return (z >= 0 && z < MESH_SIZE[2] && !point_outside_mesh(x, y) && vlmap[y][x] != NULL);
	}
	lmcell *get_lmcell(point const &p);
	void alloc(unsigned nbins, unsigned zsize, unsigned char **need_lmcell);
};


class light_source { // size = 64

	bool dynamic;
	char gl_light_id;
	float radius, radius_inv, r_inner, bwidth;
	CELL_LOC_T cent[3];
	point center;
	vector3d dir;
	colorRGBA color;

public:
	light_source() {}
	light_source(float sz, point const &p, colorRGBA const &c, bool id, vector3d const &d=plus_z, float bw=1.0, float ri=0.0, int gl_id=-1);
	void calc_cent();
	void add_color(colorRGBA const &c);
	colorRGBA const &get_color() const {return color;}
	float get_radius()           const {return radius;}
	float get_r_inner()          const {return r_inner;}
	point const &get_center()    const {return center;}
	CELL_LOC_T const *get_cent() const {return cent;}
	float get_intensity_at(point const &pos) const;
	float get_dir_intensity(vector3d const &obj_dir) const;
	bool lights_polygon(point const &pc, float rsize, vector3d const* const norm=NULL) const;
	void get_bounds(point bounds[2], int bnds[3][2], float thresh) const;
	bool is_visible()     const;
	bool is_directional() const {return (bwidth < 1.0);}
	bool is_dynamic()     const {return dynamic;}
	int  get_light_id()   const {assert(gl_light_id >= 0); return gl_light_id;}
	void shift_by(vector3d const &vd);
	void combine_with(light_source const &l);
	void draw(int ndiv) const;
	void pack_to_floatv(float *data) const;
	bool try_merge_into(light_source &ls) const;
	bool operator<(light_source const &l) const {return (radius < l.radius);} // compare radius
	bool operator>(light_source const &l) const {return (radius > l.radius);} // compare radius
};


class dls_cell {

	float z1, z2;
	vector<unsigned short> lsrc;

public:
	dls_cell() : z1(FAR_CLIP), z2(-FAR_CLIP) {}

	void clear();
	void add_light(unsigned ix, float zmin, float zmax);
	bool check_add_light(unsigned ix) const;
	size_t size() const {return lsrc.size();}
	bool empty()  const {return lsrc.empty();}
	unsigned get(unsigned i) const {return lsrc[i];} // no bounds checking
	bool check_z(float z)    const {return (!empty() && z >= z1 && z <= z2);}
	bool check_z_range(float zlo, float zhi) const {return (!empty() && zlo <= z2 && zhi >= z1);}
	void get_close_sources(point const &pos, float radius, vector<unsigned> &dlights) const;
	vector<unsigned short> const &get_src_ixs() const {return lsrc;}
};


struct flow_cache_e { // size = 16

	CELL_LOC_T f[3], t[3];
	float val;

	flow_cache_e() : val(0.0) {reset();}

	flow_cache_e(short const *f_, short const *t_, float v=0.0) : val(v) {
		f[0] = f_[0]; f[1] = f_[1]; f[2] = f_[2]; t[0] = t_[0]; t[1] = t_[1]; t[2] = t_[2];
	}
	void reset() {f[0] = f[1] = f[2] = t[0] = t[1] = t[2] = -1;}

	bool operator==(flow_cache_e const &e) const {
		return (e.f[0] == f[0] && e.f[1] == f[1] && e.f[2] == f[2] && e.t[0] == t[0] && e.t[1] == t[1] && e.t[2] == t[2]);
	}
	size_t hash() const {
		return (10831*f[0] + 12003*f[1] + 15601*f[2] + 18401*t[0] + 19813*t[1] + 21401*t[2]);
	}
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


void compute_ray_trace_lighting(unsigned ltype);


#endif

