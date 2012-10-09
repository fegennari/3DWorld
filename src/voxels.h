// 3D World - Voxel Header
// by Frank Gennari
// 2/25/12
#ifndef _VOXELS_H_
#define _VOXELS_H_

#include "3DWorld.h"
#include "model3d.h"

struct coll_tquad;


struct voxel_params_t {

	// generation parameters
	unsigned xsize, ysize, zsize, num_blocks; // num_blocks is in x and y
	float isolevel, elasticity, mag, freq, atten_thresh, tex_scale, noise_scale, noise_freq, tex_mix_saturate, z_gradient, height_eval_freq, radius_val;
	float ao_radius, ao_weight_scale, ao_atten_power;
	bool make_closed_surface, invert, remove_under_mesh, add_cobjs, normalize_to_1;
	unsigned remove_unconnected; // 0=never, 1=init only, 2=always
	unsigned atten_at_edges; // 0=no atten, 1=top only, 2=all 5 edges (excludes the bottom), 3=sphere (outer), 4=sphere (inner and outer)
	unsigned keep_at_scene_edge; // 0=don't keep, 1=always keep, 2=only when scrolling
	unsigned atten_top_mode; // 0=constant, 1=current mesh, 2=2d surface mesh
	int geom_rseed;
	
	// rendering parameters
	int texture_rseed;
	unsigned tids[2];
	colorRGBA colors[2];
	colorRGBA base_color;

	voxel_params_t() : xsize(0), ysize(0), zsize(0), num_blocks(12), isolevel(0.0), elasticity(0.5), mag(1.0), freq(1.0), atten_thresh(1.0), tex_scale(1.0), noise_scale(0.1),
		noise_freq(1.0), tex_mix_saturate(5.0), z_gradient(0.0), height_eval_freq(1.0), radius_val(0.5), ao_radius(1.0), ao_weight_scale(2.0), ao_atten_power(1.0), make_closed_surface(1),
		invert(0), remove_under_mesh(0), add_cobjs(1), normalize_to_1(1), remove_unconnected(1), atten_at_edges(0), keep_at_scene_edge(0), atten_top_mode(0), geom_rseed(123), texture_rseed(321)
	{
			tids[0] = tids[1] = 0; colors[0] = colors[1] = base_color = WHITE;
	}
};


// stored internally in yxz order
template<typename V> class voxel_grid : public vector<V> {
public:
	unsigned nx, ny, nz, xblocks, yblocks;
	vector3d vsz; // size of a voxel in x,y,z
	point center, lo_pos;

	voxel_grid() : nx(0), ny(0), nz(0), xblocks(0), yblocks(0), vsz(zero_vector) {}
	void init(unsigned nx_, unsigned ny_, unsigned nz_, vector3d const &vsz_, point const &center_, V default_val, unsigned num_blocks);
	bool is_valid_range(int i[3]) const {return (i[0] >= 0 && i[1] >= 0 && i[2] >= 0 && i[0] < (int)nx && i[1] < (int)ny && i[2] < (int)nz);}

	void get_xyz(point const &p, int xyz[3]) const { // returns whether or not the point was inside the voxel volume
		UNROLL_3X(xyz[i_] = int((p[i_] - lo_pos[i_])/vsz[i_]);); // convert to voxel space
	}
	bool get_ix(point const &p, unsigned &ix) const { // returns whether or not the point was inside the voxel volume
		int i[3]; // x,y,z
		get_xyz(p, i);
		if (!is_valid_range(i)) return 0;
		ix = i[2] + (i[0] + i[1]*nx)*nz;
		return 1;
	}
	unsigned get_ix(unsigned x, unsigned y, unsigned z) const {
		//assert(x < nx && y < ny && z < nz);
		return (z + (x + y*nx)*nz);
	}
	point get_pt_at(unsigned x, unsigned y, unsigned z) const  {return (point(x, y, z)*vsz + lo_pos);}
	V const &get   (unsigned x, unsigned y, unsigned z) const  {return operator[](get_ix(x, y, z));}
	V &get_ref     (unsigned x, unsigned y, unsigned z)        {return operator[](get_ix(x, y, z));}
	void set       (unsigned x, unsigned y, unsigned z, V const &val) {operator[](get_ix(x, y, z)) = val;}
	cube_t get_raw_bbox() const {return cube_t(lo_pos, center + (center - lo_pos));}
};

typedef voxel_grid<float> float_voxel_grid;


class voxel_manager : public float_voxel_grid {

protected:
	bool use_mesh;
	voxel_params_t params;
	voxel_grid<unsigned char> outside;

	point interpolate_pt(float isolevel, point const &pt1, point const &pt2, float const val1, float const val2) const;
	void calc_outside_val(unsigned x, unsigned y, unsigned z, bool is_under_mesh);
	void remove_unconnected_outside_range(bool keep_at_edge, unsigned x1, unsigned y1, unsigned x2, unsigned y2,
		vector<unsigned> *xy_updated, vector<point> *updated_pts);
	unsigned get_triangles_for_voxel(vector<triangle> &triangles, unsigned x, unsigned y, unsigned z, bool count_only) const;
	void add_cobj_voxels(coll_obj &cobj, float filled_val);

public:
	voxel_manager(bool use_mesh_=0) : use_mesh(use_mesh_) {}
	void set_params(voxel_params_t const &p) {params = p;}
	void clear();
	void create_procedural(float mag, float freq, vector3d const &offset, bool normalize_to_1, int rseed1, int rseed2);
	void create_from_cobjs(coll_obj_group &cobjs, float filled_val=1.0);
	void atten_at_edges(float val);
	void atten_at_top_only(float val);
	void atten_to_sphere(float val, float inner_radius, bool atten_inner);
	void determine_voxels_outside();
	void remove_unconnected_outside();
	bool is_outside(unsigned ix) const {assert(ix < outside.size()); return((outside[ix]&3) != 0);}
	bool point_inside_volume(point const &pos) const;
	bool point_intersect(point const &center, point *int_pt) const;
	bool sphere_intersect(point const &center, float radius, point *int_pt) const;
	bool line_intersect(point const &p1, point const &p2, point *int_pt) const;
	unsigned upload_to_3d_texture(int wrap) const;
};


class noise_texture_manager_t {

	unsigned noise_tid, tsize;
	voxel_manager voxels;

public:
	noise_texture_manager_t() : noise_tid(0), tsize(0) {}
	void setup(unsigned size, int rseed=321, float mag=1.0, float freq=1.0, vector3d const &offset=zero_vector);
	void bind_texture(unsigned tu_id) const;
	void clear();
	float eval_at(point const &pos) const;
};


class voxel_model : public voxel_manager {

protected:
	bool volume_added;
	typedef vert_norm vertex_type_t;
	typedef vntc_vect_block_t<vertex_type_t> tri_data_t;
	tri_data_t tri_data;
	noise_texture_manager_t noise_tex_gen;
	std::set<unsigned> modified_blocks;
	voxel_grid<unsigned char> ao_lighting;

	struct step_dir_t {
		unsigned nsteps;
		int dir[3], dist_per_step;
		step_dir_t(int x, int y, int z, unsigned n, int d) : nsteps(n), dist_per_step(d) {dir[0] = x; dir[1] = y; dir[2] = z;}
	};
	vector<step_dir_t> ao_dirs;

	struct pt_ix_t {
		point pt;
		unsigned ix;
		pt_ix_t(point const &pt_=all_zeros, unsigned const ix_=0) : pt(pt_), ix(ix_) {}
	};

	vector<pt_ix_t> pt_to_ix;

	struct comp_by_dist {
		point const p;
		comp_by_dist(point const &p_) : p(p_) {}
		bool operator()(pt_ix_t const &a, pt_ix_t const &b) const {return (p2p_dist_sq(a.pt, p) < p2p_dist_sq(b.pt, p));}
	};

	void remove_unconnected_outside_modified_blocks(void);
	unsigned get_block_ix(unsigned voxel_ix) const;
	virtual bool clear_block(unsigned block_ix);
	unsigned create_block(unsigned block_ix, bool first_create, bool count_only);
	void calc_ao_dirs();
	virtual void calc_ao_lighting_for_block(unsigned block_ix, bool increase_only);
	void calc_ao_lighting();

	virtual void maybe_create_fragments(point const &center, float radius, int shooter, unsigned num_fragments) const {} // do nothing
	virtual void create_block_hook(unsigned block_ix, vector<triangle> const &triangles) {}
	virtual void update_blocks_hook(vector<unsigned> const &blocks_to_update, unsigned num_added) {}
	virtual void pre_build_hook() {}

public:
	voxel_model(bool use_mesh_) : voxel_manager(use_mesh_), volume_added(0) {}
	virtual ~voxel_model() {}
	void clear();
	bool update_voxel_sphere_region(point const &center, float radius, float val_at_center, point *damage_pos=NULL, int shooter=-1, unsigned num_fragments=0);
	unsigned get_texture_at(point const &pos) const;
	void proc_pending_updates();
	void build(bool verbose);
	virtual void setup_tex_gen_for_rendering(shader_t &s);
	void core_render(shader_t &s, bool is_shadow_pass);
	void render(bool is_shadow_pass);
	virtual void free_context();
	float eval_noise_texture_at(point const &pos) const;
	float get_ao_lighting_val(point const &pos) const;
	cube_t get_bcube() const {return ((tri_data.empty()) ? cube_t(center, center) : tri_data.get_bbox());}
	sphere_t get_bsphere() const;
};


class voxel_model_ground : public voxel_model {

	bool add_cobjs, add_as_fixed;
	vector<unsigned> last_blocks_updated;

	struct data_block_t {
		vector<int> cids; // references into coll_objects
		//unsigned tri_data_ix;
		void clear() {cids.clear();}
	};
	vector<data_block_t> data_blocks;

	virtual bool clear_block(unsigned block_ix);
	virtual void maybe_create_fragments(point const &center, float radius, int shooter, unsigned num_fragments) const;
	virtual void create_block_hook(unsigned block_ix, vector<triangle> const &triangles);
	virtual void update_blocks_hook(vector<unsigned> const &blocks_to_update, unsigned num_added);
	virtual void pre_build_hook();

public:
	voxel_model_ground() : voxel_model(1), add_cobjs(0), add_as_fixed(0) {}
	void clear();
	void build(bool add_cobjs_, bool add_as_fixed_, bool verbose);
};


class voxel_model_space : public voxel_model {

	unsigned ao_tid, shadow_tid;
	vector<triangle> shadow_edge_tris;

	void free_ao_and_shadow_texture() {free_texture(ao_tid); free_texture(shadow_tid);}
	virtual void calc_ao_lighting_for_block(unsigned block_ix, bool increase_only);
	void calc_shadows(voxel_grid<unsigned char> &shadow_data) const;
	void extract_shadow_edges(voxel_grid<unsigned char> const &shadow_data);

public:
	voxel_model_space() : voxel_model(0), ao_tid(0), shadow_tid(0) {}
	void clear() {voxel_model::clear(); shadow_edge_tris.clear();}
	virtual void free_context() {voxel_model::free_context(); free_ao_and_shadow_texture();}
	virtual void setup_tex_gen_for_rendering(shader_t &s);
	vector<triangle> const &get_shadow_edge_tris() const {return shadow_edge_tris;}
};


#endif

