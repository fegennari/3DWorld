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

#define ADD_LIGHT_CONTRIB(c, C) C[0] += c[0]; C[1] += c[1]; C[2] += c[2];


struct lmcell { // size = 40

	float c[3], ac[3], v, smoke; // c: RGB, ac: indirect lighting color RGB
	unsigned char lflow[3], pflow[3]; // flow: x, y, z
	
	lmcell() : v(0.0), smoke(0.0) {
		UNROLL_3X(c[i_] = ac[i_] = 0.0; lflow[i_] = pflow[i_] = 255;)
	}
	// c[0],c[1],c[2] : ac[0],ac[1],ac[2],v
	float       *get_offset(bool local)       {return (local ? c : ac);}
	float const *get_offset(bool local) const {return (local ? c : ac);}
	unsigned     get_dsz   (bool local) const {return (local ? 3 : 4 );}
};


class lmap_manager_t {

	vector<lmcell> vldata_alloc;
	unsigned lm_zsize;

public:
	lmcell ***vlmap; // y, x, z

	lmap_manager_t() : lm_zsize(0), vlmap(NULL) {}
	void clear() {vldata_alloc.clear();} // reset vlmap to NULL?
	size_t size() const {return vldata_alloc.size();}
	bool read_data_from_file(char const *const fn, bool local);
	bool write_data_to_file(char const *const fn, bool local) const;
	void global_light_scale(float scale);
	void local_light_scale(float scale);

	inline bool is_valid_cell(int x, int y, int z) const {
		return (z >= 0 && z < MESH_SIZE[2] && !point_outside_mesh(x, y) && vlmap[y][x] != NULL);
	}
	lmcell *get_lmcell(point const &p);
	void alloc(unsigned nbins, unsigned zsize, unsigned char **need_lmcell);
	void normalize_light_val(float min_light, float max_light, float light_scale, float light_off);
};


class light_source { // size = 64

	bool dynamic;
	float radius, radius_inv, r_inner, bwidth;
	CELL_LOC_T cent[3];
	point center;
	vector3d dir;
	colorRGBA color;

public:
	light_source() {}
	light_source(float sz, point const &p, colorRGBA const &c, bool id, vector3d const &d=plus_z, float bw=1.0, float ri=0.0);
	void calc_cent();
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
	void shift_by(vector3d const &vd);
	void combine_with(light_source const &l);
	void draw(int ndiv) const;
	void pack_to_floatv(float *data) const;
	
	bool operator<(light_source const &l) { // compare: y, x, z, r
		if (center.y < l.center.y) return 1;
		if (center.y > l.center.y) return 0;
		if (center.x < l.center.x) return 1;
		if (center.x > l.center.x) return 0;
		if (center.z < l.center.z) return 1;
		if (center.z > l.center.z) return 0;
		return (radius < l.radius);
	}
};


class dls_cell {

	float z1, z2;
	vector<unsigned> lsrc;

public:
	dls_cell() : z1(FAR_CLIP), z2(-FAR_CLIP) {}

	void clear() {
		if (lsrc.capacity() > INIT_CCELL_SIZE) lsrc.clear(); else lsrc.resize(0);
		z1 =  FAR_CLIP;
		z2 = -FAR_CLIP;
	}
	void add_light(unsigned ix, float zmin, float zmax) {
		if (lsrc.capacity() == 0) lsrc.reserve(INIT_CCELL_SIZE);
		lsrc.push_back(ix);
		z1 = min(z1, zmin);
		z2 = max(z2, zmax);
	}
	size_t size() const {return lsrc.size();}
	bool empty()  const {return lsrc.empty();}
	unsigned get(unsigned i) const {return lsrc[i];} // no bounds checking
	bool check_z(float z)    const {return (!empty() && z >= z1 && z <= z2);}
	bool check_z_range(float zlo, float zhi) const {return (!empty() && zlo <= z2 && zhi >= z1);}
	void get_close_sources(point const &pos, float radius, vector<unsigned> &dlights) const;
	vector<unsigned> const &get_src_ixs() const {return lsrc;}
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


bool is_under_mesh(point const &p);
void compute_ray_trace_lighting_global();
void compute_ray_trace_lighting_local();


#endif

