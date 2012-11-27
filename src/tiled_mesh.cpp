// 3D World - Tiled Landscape Mesh Generation/Drawing Code
// by Frank Gennari
// 9/26/10

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "small_tree.h"
#include "scenery.h"
#include "grass.h"
#include "tree_3dw.h"


bool const DEBUG_TILES        = 0;
bool const DEBUG_TILE_BOUNDS  = 0;
bool const ENABLE_TREE_LOD    = 1; // faster but has popping artifacts
bool const ENABLE_TERRAIN_ENV = 1;
bool const GRASS_CLOUD_SHADOWS= 1; // slow, but looks nice
int  const TILE_RADIUS        = 6; // in mesh sizes
unsigned const NUM_LODS       = 5; // > 0
unsigned const NORM_TEXELS    = 512;
float const FOG_DIST_TILES    = 1.4;
float const DRAW_DIST_TILES   = 1.45;
float const CREATE_DIST_TILES = 1.5;
float const CLEAR_DIST_TILES  = 1.5;
float const DELETE_DIST_TILES = 1.7;
float const TREE_LOD_THRESH   = 5.0;
float const GEOMORPH_THRESH   = 5.0;
float const SCENERY_THRESH    = 1.6;
float const GRASS_THRESH      = 1.5;
float const GRASS_LOD_SCALE   = 16.0;


extern bool inf_terrain_scenery;
extern unsigned grass_density;
extern int xoff, yoff, island, DISABLE_WATER, display_mode, show_fog, tree_mode, leaf_color_changed, ground_effects_level;
extern float zmax, zmin, water_plane_z, mesh_scale, mesh_scale_z, vegetation, relh_adj_tex, grass_length, grass_width;
extern point sun_pos, moon_pos;
extern float h_dirt[];
extern texture_t textures[];
extern tree_data_manager_t tree_data_manager;

bool enable_terrain_env(ENABLE_TERRAIN_ENV);


float get_scaled_tile_radius  () {return TILE_RADIUS*(X_SCENE_SIZE + Y_SCENE_SIZE);}
float get_inf_terrain_fog_dist() {return FOG_DIST_TILES*get_scaled_tile_radius();}
bool is_water_enabled     () {return (!DISABLE_WATER && (display_mode & 0x04) != 0);}
bool pine_trees_enabled   () {return ((tree_mode & 2) && vegetation > 0.0);}
bool decid_trees_enabled  () {return ((tree_mode & 1) && vegetation > 0.0);}
bool any_trees_enabled    () {return (pine_trees_enabled() || decid_trees_enabled());}
bool scenery_enabled      () {return (inf_terrain_scenery && SCENERY_THRESH > 0.0);}
bool is_grass_enabled     () {return ((display_mode & 0x02) && GRASS_THRESH > 0.0 && grass_density > 0);}
bool cloud_shadows_enabled() {return (ground_effects_level >= 2);}
float get_tiled_terrain_water_level() {return (is_water_enabled() ? water_plane_z : zmin);}


grass_tile_manager_t grass_tile_manager;

void update_tiled_terrain_grass_vbos() {
	grass_tile_manager.invalidate_vbo();
}


void bind_texture_tu(unsigned tid, unsigned tu_id) {

	assert(tid);
	set_multitex(tu_id);
	bind_2d_texture(tid);
	set_multitex(0);
}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator<(tile_xy_pair const &t) const {return ((y == t.y) ? (x < t.x) : (y < t.y));}
};


#define BILINEAR_INTERP(arr, var, x, y) (y*(x*arr[1][1].var + (1.0-x)*arr[1][0].var) + (1.0-y)*(x*arr[0][1].var + (1.0-x)*arr[0][0].var))


class tile_t;
tile_t *get_tile_from_xy(tile_xy_pair const &tp);


class tile_t {

	struct offset_t {
		int dxoff, dyoff;

		offset_t(int dxoff_=0, int dyoff_=0) : dxoff(dxoff_), dyoff(dyoff_) {}
		void set_from_xyoff2() {dxoff = -xoff2; dyoff = -yoff2;}
		vector3d get_xlate() const {return vector3d(((xoff - xoff2) - dxoff)*DX_VAL, ((yoff - yoff2) - dyoff)*DY_VAL, 0.0);}
		vector3d subtract_from(offset_t const &o) const {return vector3d((o.dxoff - dxoff)*DX_VAL, (o.dyoff - dyoff)*DY_VAL, 0.0);}
	};

	int x1, y1, x2, y2;
	unsigned weight_tid, height_tid, shadow_normal_tid, vbo, ivbo[NUM_LODS];
	unsigned size, stride, zvsize, base_tsize, gen_tsize;
	float radius, mzmin, mzmax, ptzmax, dtzmax, trmax, xstart, ystart, xstep, ystep;
	bool shadows_invalid, weights_invalid, in_queue;
	offset_t mesh_off, ptree_off, dtree_off, scenery_off;
	vector<float> zvals;
	vector<unsigned char> tree_map;
	vector<unsigned char> smask[NUM_LIGHT_SRC];
	vector<float> sh_out[NUM_LIGHT_SRC][2];
	small_tree_group pine_trees;
	scenery_group scenery;
	tree_cont_t decid_trees;

	struct grass_block_t {
		unsigned ix; // 0 is unused
		float zmin, zmax;
		grass_block_t() : ix(0), zmin(0.0), zmax(0.0) {}
	};

	vector<grass_block_t> grass_blocks;

	struct terrain_params_t {
		float hoff, hscale, veg, grass, dirt;
		terrain_params_t() : hoff(0.0), hscale(1.0), veg(1.0), grass(1.0), dirt(0.0) {}
	};

	terrain_params_t params[2][2]; // {ylo,yhi} x {xlo,xhi}

	void update_terrain_params() {
		float const off_mult(0.4), height_mult(0.8), dirt_mult(1.0), veg_mult(5.0), off_scale(1.0);

		for (unsigned yp = 0; yp < 2; ++yp) {
			for (unsigned xp = 0; xp < 2; ++xp) {
				terrain_params_t &param(params[yp][xp]);
				float const xv(mesh_scale*get_xval(xp ? x2 : x1)), yv(mesh_scale*get_yval(yp ? y2 : y1));
				//param.hoff   = off_scale*eval_mesh_sin_terms(off_mult*xv+123, off_mult*yv+456);
				//param.hscale = min(2.0f, max(0.5f, 0.5f*fabs(eval_mesh_sin_terms(height_mult*xv+789, height_mult*yv+111))));
				float const veg_val(eval_mesh_sin_terms(veg_mult*xv, veg_mult*yv));
				param.veg    = CLIP_TO_01(5.000f*(veg_val + 1.5f));
				param.grass  = CLIP_TO_01(100.0f*(veg_val + 3.0f)); // depends on hoff?
				param.dirt   = CLIP_TO_01(5.0f*(eval_mesh_sin_terms(dirt_mult*xv, dirt_mult*yv) + 1.0f));
			}
		}
	}

public:
	typedef point vert_type_t;

	tile_t() : weight_tid(0), height_tid(0), shadow_normal_tid(0), vbo(0), size(0), stride(0), zvsize(0), gen_tsize(0), decid_trees(tree_data_manager) {
		init_vbo_ids();
	}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~tile_t() {clear_vbo_tid(1,1);}
	
	tile_t(unsigned size_, int x, int y) : weight_tid(0), height_tid(0), shadow_normal_tid(0), vbo(0), size(size_), stride(size+1), zvsize(stride+1),
		gen_tsize(0), trmax(0.0), shadows_invalid(1), weights_invalid(1), in_queue(0), mesh_off(xoff-xoff2, yoff-yoff2), decid_trees(tree_data_manager)
	{
		assert(size > 0);
		x1 = x*size;
		y1 = y*size;
		x2 = x1 + size;
		y2 = y1 + size;
		calc_start_step(0, 0);
		radius = calc_radius();
		mzmin  = mzmax = ptzmax = dtzmax = get_camera_pos().z;
		base_tsize = NORM_TEXELS;
		init_vbo_ids();
	}
	void invalidate_shadows() {shadows_invalid = 1;}
	float calc_radius() const {return 0.5*sqrt(xstep*xstep + ystep*ystep)*size;} // approximate (lower bound)
	float get_zmin() const {return mzmin;}
	float get_zmax() const {return mzmax;}
	float get_tile_zmax() const {return max((mzmax + (grass_blocks.empty() ? 0.0f : grass_length)), max(ptzmax, dtzmax));}
	bool has_water() const {return (mzmin < water_plane_z);}
	bool all_water() const {return (mzmax < water_plane_z);} // get_tile_zmax()?
	bool pine_trees_generated() const {return pine_trees.generated;}
	bool has_pine_trees() const {return (pine_trees_generated() && !pine_trees.empty());}
	float get_avg_veg() const {return 0.25*(params[0][0].veg + params[0][1].veg + params[1][0].veg + params[1][1].veg);}

	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5*(mzmin + mzmax));
	}
	cube_t get_cube() const {
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2)), z2(get_tile_zmax());
		return cube_t(xv1-trmax, xv1+(x2-x1)*DX_VAL+trmax, yv1-trmax, yv1+(y2-y1)*DY_VAL+trmax, mzmin, z2);
	}
	float get_min_dist_to_pt(point const &pt) const {
		point p1(pt), p2(get_center());
		bool const ret(do_line_clip(p1, p2, get_cube().d)); // only clip in x and y?
		assert(ret);
		return p2p_dist(pt, p1);
	}
	float get_max_xy_dist_to_pt(point const &pt) const {
		cube_t const bc(get_cube());
		float const dx(max(fabs(pt.x - bc.d[0][0]), fabs(pt.x - bc.d[0][1])));
		float const dy(max(fabs(pt.y - bc.d[1][0]), fabs(pt.y - bc.d[1][1])));
		return sqrt(dx*dx + dy*dy);
	}
	bool contains_camera() const {
		return get_cube().contains_pt_xy(get_camera_pos());
	}

	void calc_start_step(int dx, int dy) {
		xstart = get_xval(x1 + dx);
		ystart = get_yval(y1 + dy);
		xstep  = (get_xval(x2 + dx) - xstart)/size;
		ystep  = (get_yval(y2 + dy) - ystart)/size;
	}

	unsigned get_gpu_mem() const {
		unsigned mem(pine_trees.get_gpu_mem() + decid_trees.get_gpu_mem() + scenery.get_gpu_mem());
		unsigned const num_texels(stride*stride);
		if (vbo > 0) mem += 2*stride*size*sizeof(vert_type_t);
		if (weight_tid > 0) mem += 4*num_texels; // 4 bytes per texel (RGBA8)
		if (height_tid > 0) mem += 2*num_texels; // 2 bytes per texel (L16)
		if (shadow_normal_tid > 0) mem += 4*num_texels; // 4 bytes per texel (L8)

		for (unsigned i = 0; i < NUM_LODS; ++i) {
			if (ivbo[i] > 0) mem += (size>>i)*(size>>i)*sizeof(unsigned short);
		}
		return mem;
	}
	unsigned get_tree_mem() const { // only accounts for top-level class memory
		return (pine_trees.capacity()*sizeof(small_tree) + decid_trees.capacity()*sizeof(tree));
	}

	void clear() {
		clear_vbo_tid(1, 1);
		tree_map.clear();
		zvals.clear();
		clear_shadows();
		pine_trees.clear_all();
		decid_trees.clear();
		scenery.clear();
		grass_blocks.clear();
	}

	void clear_shadows() {
		for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
			smask[l].clear();
			for (unsigned d = 0; d < 2; ++d) sh_out[l][d].clear();
		}
		invalidate_shadows();
	}

	void clear_tids() {
		free_texture(weight_tid);
		free_texture(height_tid);
		free_texture(shadow_normal_tid);
		gen_tsize = 0;
	}

	void init_vbo_ids() {
		vbo = 0;
		for (unsigned i = 0; i < NUM_LODS; ++i) {ivbo[i] = 0;}
	}

	void clear_vbo_tid(bool vclear, bool tclear) {
		if (vclear) {
			delete_vbo(vbo);
			for (unsigned i = 0; i < NUM_LODS; ++i) {delete_vbo(ivbo[i]);}
			init_vbo_ids();
			clear_shadows();
			pine_trees.clear_vbos();
			decid_trees.clear_vbos(); // only necessary if not using instancing
			scenery.clear_vbos_and_dlists();
		}
		if (tclear) {clear_tids();}
	}

	void create_xy_arrays(mesh_xy_grid_cache_t &height_gen, unsigned xy_size, float xy_scale) {
		height_gen.build_arrays((xy_scale*(xstart - 0.5*xstep)), (xy_scale*(ystart + 0.5*ystep)), xy_scale*xstep, xy_scale*ystep, xy_size, xy_size);
	}

	void create_zvals(mesh_xy_grid_cache_t &height_gen) {
		//RESET_TIME;
		if (enable_terrain_env) {update_terrain_params();}
		zvals.resize(zvsize*zvsize);
		calc_start_step(0, 0);
		create_xy_arrays(height_gen, zvsize, 1.0);
		mzmin =  FAR_CLIP;
		mzmax = -FAR_CLIP;
		float const xy_mult(1.0/float(size));

		#pragma omp parallel for schedule(static,1)
		for (int y = 0; y < (int)zvsize; ++y) {
			for (unsigned x = 0; x < zvsize; ++x) {
				float const xv(float(x)*xy_mult), yv(float(y)*xy_mult);
				float const hoff(BILINEAR_INTERP(params, hoff, xv, yv)), hscale(BILINEAR_INTERP(params, hscale, xv, yv));
				zvals[y*zvsize + x] = hoff + hscale*height_gen.eval_index(x, y);
			}
		}
		for (vector<float>::const_iterator i = zvals.begin(); i != zvals.end(); ++i) {
			mzmin = min(mzmin, *i);
			mzmax = max(mzmax, *i);
		}
		assert(mzmin <= mzmax);
		radius = 0.5*sqrt((xstep*xstep + ystep*ystep)*size*size + (mzmax - mzmin)*(mzmax - mzmin));
		//PRINT_TIME("Create Zvals");
		if (DEBUG_TILES) {cout << "new tile coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;}
	}

	inline vector3d get_norm(unsigned ix) const {
		return vector3d(DY_VAL*(zvals[ix] - zvals[ix + 1]), DX_VAL*(zvals[ix] - zvals[ix + zvsize]), dxdy).get_norm();
	}

	// *** shadows ***

	void calc_shadows_for_light(unsigned l) {
		assert(!smask[l].empty());
		tile_xy_pair const tp(x1/(int)size, y1/(int)size);

		// pull from adjacent tiles that already had their shadows calculated
		float const *sh_in[2] = {0, 0};
		point const lpos(get_light_pos(l));
		tile_xy_pair const adj_tp[2] = {tile_xy_pair((tp.x + ((lpos.x < 0.0) ? -1 : 1)), tp.y),
										tile_xy_pair(tp.x, (tp.y + ((lpos.y < 0.0) ? -1 : 1)))}; // toward the light source

		for (unsigned d = 0; d < 2; ++d) { // d = tile adjacency dimension, shared edge is in !d
			sh_out[l][!d].resize(zvsize, MESH_MIN_Z); // init value really should not be used, but it sometimes is
			tile_t *adj_tile(get_tile_from_xy(adj_tp[d]));
			if (adj_tile == NULL) continue; // no adjacent tile
			vector<float> const &adj_sh_out(adj_tile->sh_out[l][!d]);
				
			if (adj_sh_out.empty()) { // adjacent tile not initialized
				adj_tile->calc_shadows((l == LIGHT_SUN), (l == LIGHT_MOON), 1); // recursive call on adjacent tile
			}
			assert(adj_sh_out.size() == zvsize);
			sh_in[!d] = &adj_sh_out.front(); // chain our input to our neighbor's output
		}

		// calculate shadows of current tile
		calc_mesh_shadows(l, lpos, &zvals.front(), &smask[l].front(), zvsize, zvsize,
			sh_in[0], sh_in[1], &sh_out[l][0].front(), &sh_out[l][1].front());
		invalidate_shadows();
	}

	static void proc_tile_queue(tile_t *init_tile, unsigned l) {
		point const lpos(get_light_pos(l));
		deque<tile_t *> tile_queue;
		tile_queue.push_front(init_tile);
		init_tile->in_queue = 1;

		while (!tile_queue.empty()) {
			tile_t *t(tile_queue.back());
			tile_queue.pop_back();
			assert(t->in_queue);
			t->in_queue = 0;
			vector<float> const prev_sh_out[2] = {t->sh_out[l][0], t->sh_out[l][1]};
			t->calc_shadows_for_light(l);
			tile_xy_pair const tp(t->x1/int(t->size), t->y1/int(t->size));
			tile_xy_pair const adj_tp2[2] = {tile_xy_pair((tp.x + ((lpos.x < 0.0) ? 1 : -1)), tp.y),
											 tile_xy_pair(tp.x, (tp.y + ((lpos.y < 0.0) ? 1 : -1)))}; // away from the light source

			for (unsigned d = 0; d < 2; ++d) { // d = tile adjacency dimension, shared edge is in !d
				if (t->sh_out[l][!d] == prev_sh_out[!d]) continue; // unchanged, no update needed
				tile_t *adj_tile(get_tile_from_xy(adj_tp2[d]));
				if (adj_tile == NULL || adj_tile->smask[l].empty() || adj_tile->in_queue) continue; // no adjacent tile, not initialized, or already in queue
				tile_queue.push_front(adj_tile); // changed, push to adjacent tiles
				adj_tile->in_queue = 1;
			}
		}
	}

	void calc_shadows(bool calc_sun, bool calc_moon, bool no_push=0) {
		bool calc_light[NUM_LIGHT_SRC] = {0};
		calc_light[LIGHT_SUN ] = calc_sun;
		calc_light[LIGHT_MOON] = calc_moon;

		for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) { // calculate mesh shadows for each light source
			if (!calc_light[l])    continue; // light not enabled
			if (!smask[l].empty()) continue; // already calculated (cached)
			smask[l].resize(zvals.size());

			if (no_push) {
				calc_shadows_for_light(l);
			}
			else {
				proc_tile_queue(this, l);
			}
		}
	}

	tile_t *get_adj_tile_smap(int dx, int dy) const {
		tile_xy_pair const tp((x1/(int)size)+dx, (y1/(int)size)+dy);
		tile_t *adj_tile(get_tile_from_xy(tp));
		return ((adj_tile && !adj_tile->tree_map.empty()) ? adj_tile : NULL);
	}

	void push_tree_ao_shadow(int dx, int dy, point const &pos, float radius) const {
		tile_t *const adj_tile(get_adj_tile_smap(dx, dy));
		if (!adj_tile) return;
		point const pos2(pos + mesh_off.subtract_from(adj_tile->mesh_off));
		adj_tile->add_tree_ao_shadow(pos2, radius, 1);
	}

	bool add_tree_ao_shadow(point const &pos, float radius, bool no_adj_test) {
		int const xc(round_fp((pos.x - xstart)/xstep)), yc(round_fp((pos.y - ystart)/ystep));
		int rval(max(int(radius/xstep), int(radius/ystep)) + 1);
		int const x1(max(0, xc-rval)), y1(max(0, yc-rval)), x2(min((int)size, xc+rval)), y2(min((int)size, yc+rval));
		bool on_edge(0);
		float const scale(0.6/rval);

		for (int y = y1; y <= y2; ++y) {
			for (int x = x1; x <= x2; ++x) {
				float const dx(abs(x - xc)), dy(abs(y - yc)), dist(sqrt(dx*dx + dy*dy));
				if (dist < rval) {tree_map[y*stride + x] *= (0.2 + 0.8*scale*dist);}
			}
		}
		if (!no_adj_test) {
			bool const x_test[3] = {(xc <= rval), 1, (xc >= (int)size-rval)};
			bool const y_test[3] = {(yc <= rval), 1, (yc >= (int)size-rval)};

			for (int dy = -1; dy <= 1; ++dy) {
				for (int dx = -1; dx <= 1; ++dx) {
					if (dx == 0 && dy == 0) continue;
					if (x_test[dx+1] && y_test[dy+1]) {push_tree_ao_shadow(dx, dy, pos, radius); on_edge = 1;}
				}
			}
		}
		invalidate_shadows();
		weights_invalid = 1; // Note: may be slow, and doesn't have a big impact
		return on_edge;
	}

	void apply_ao_shadows_for_trees(tile_t const *const tile, bool no_adj_test) {
		if (pine_trees_enabled()) {
			assert(pine_trees_generated());
			point const pt_off(tile->ptree_off.subtract_from(mesh_off));

			for (small_tree_group::const_iterator i = tile->pine_trees.begin(); i != tile->pine_trees.end(); ++i) {
				add_tree_ao_shadow((i->get_pos() + pt_off), i->get_pine_tree_radius(), no_adj_test);
			}
		}
		if (decid_trees_enabled()) {
			assert(decid_trees.was_generated());
			point const dt_off(tile->dtree_off.subtract_from(mesh_off));

			for (tree_cont_t::const_iterator i = tile->decid_trees.begin(); i != tile->decid_trees.end(); ++i) {
				add_tree_ao_shadow((i->get_center() + dt_off), 0.6*i->get_radius(), no_adj_test); // less dense => smaller radius
			}
		}
		if (!no_adj_test) { // pull mode
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dx = -1; dx <= 1; ++dx) {
					if (dx == 0 && dy == 0) continue;
					tile_t const *const adj_tile(get_adj_tile_smap(dx, dy));
					if (adj_tile) {apply_ao_shadows_for_trees(adj_tile, 1);}
				}
			}
		}
	}

	void apply_tree_ao_shadows() { // should this generate a float or unsigned char shadow weight instead?
		tree_map.resize(stride*stride, 255);
		apply_ao_shadows_for_trees(this, 0);
	}

	struct norm_comp_with_shadow {
		unsigned char v[4];
	};

	void check_shadow_map_and_normal_texture() {
		if (shadows_invalid) {free_texture(shadow_normal_tid);}
		if (shadow_normal_tid) return; // up-to-date
		setup_texture(shadow_normal_tid, GL_MODULATE, 0, 0, 0, 0, 0);
		vector<norm_comp_with_shadow> data(stride*stride);
		bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6);
		assert(has_sun || has_moon);
		//RESET_TIME;
		calc_shadows(has_sun, has_moon);
		//PRINT_TIME("Calc Shadows");

		for (unsigned y = 0; y < stride; ++y) {
			for (unsigned x = 0; x < stride; ++x) {
				unsigned const ix(y*stride + x), ix2(y*zvsize + x);
				vector3d const norm(get_norm(y*zvsize + x));
				UNROLL_3X(data[ix].v[i_] = (unsigned char)(127.0*(norm[i_] + 1.0)););
				unsigned char shadow_val(tree_map.empty() ? 255 : tree_map[ix]);

				if (has_sun && has_moon) {
					bool const no_sun( (smask[LIGHT_SUN ][ix2] & SHADOWED_ALL) != 0);
					bool const no_moon((smask[LIGHT_MOON][ix2] & SHADOWED_ALL) != 0);
					shadow_val *= blend_light(light_factor, !no_sun, !no_moon);
				}
				else if (smask[has_sun ? LIGHT_SUN : LIGHT_MOON][ix2] & SHADOWED_ALL) {
					shadow_val = 0; // full shadow
				}
				data[ix].v[3] = shadow_val;
			}
		}
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, stride, stride, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
		shadows_invalid = 0;
	}

	// *** mesh creation ***

	void create_data(vector<vert_type_t> &data, vector<unsigned short> indices[NUM_LODS]) {
		//RESET_TIME;
		assert(zvals.size() == zvsize*zvsize);
		calc_start_step(mesh_off.dxoff, mesh_off.dyoff);
		data.resize(stride*stride);

		for (unsigned y = 0; y <= size; ++y) {
			for (unsigned x = 0; x <= size; ++x) {
				data[y*stride + x].assign((xstart + x*xstep), (ystart + y*ystep), zvals[y*zvsize + x]);
			}
		}
		for (unsigned i = 0; i < NUM_LODS; ++i) {
			indices[i].resize(4*(size>>i)*(size>>i));
			unsigned const step(1 << i);

			for (unsigned y = 0; y < size; y += step) {
				for (unsigned x = 0; x < size; x += step) {
					unsigned const vix(y*stride + x), iix(4*((y>>i)*(size>>i) + (x>>i)));
					indices[i][iix+0] = vix;
					indices[i][iix+1] = vix + step*stride;
					indices[i][iix+2] = vix + step*stride + step;
					indices[i][iix+3] = vix + step;
				}
			}
		}
		//PRINT_TIME("Create Data");
	}

	unsigned get_grass_block_dim() const {return (1+(size-1)/GRASS_BLOCK_SZ);} // ceil

	void create_texture(mesh_xy_grid_cache_t &height_gen) {
		assert(weight_tid == 0);
		assert(!island);
		assert(zvals.size() == zvsize*zvsize);
		//RESET_TIME;
		weights_invalid = 0;
		grass_blocks.clear();
		unsigned const grass_block_dim(get_grass_block_dim()), tsize(stride);
		unsigned char *data(new unsigned char[4*tsize*tsize]); // RGBA
		float const xy_mult(1.0/float(size)), water_level(get_water_z_height());
		float const MESH_NOISE_SCALE = 0.003;
		float const MESH_NOISE_FREQ  = 80.0;
		float const dz_inv(1.0/(zmax - zmin)), noise_scale(MESH_NOISE_SCALE*mesh_scale_z);
		int k1, k2, k3, k4, sand_tex_ix(-1), dirt_tex_ix(-1), grass_tex_ix(-1), rock_tex_ix(-1), snow_tex_ix(-1);
		float t;

		for (unsigned i = 0; i < NTEX_DIRT; ++i) {
			switch (lttex_dirt[i].id) {
			case SAND_TEX:   sand_tex_ix  = i; break;
			case DIRT_TEX:   dirt_tex_ix  = i; break;
			case GROUND_TEX: grass_tex_ix = i; break;
			case ROCK_TEX:   rock_tex_ix  = i; break;
			case SNOW_TEX:   snow_tex_ix  = i; break;
			default: assert(0);
			}
		}
		assert(sand_tex_ix >= 0 && dirt_tex_ix >= 0 && grass_tex_ix >= 0 && rock_tex_ix >= 0 && snow_tex_ix >= 0);
		create_xy_arrays(height_gen, zvsize, MESH_NOISE_FREQ);

		//#pragma omp parallel for schedule(static,1)
		for (unsigned y = 0; y < tsize-DEBUG_TILE_BOUNDS; ++y) {
			for (unsigned x = 0; x < tsize-DEBUG_TILE_BOUNDS; ++x) {
				float weights[NTEX_DIRT] = {0};
				unsigned const ix(y*zvsize + x);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const rand_offset(noise_scale*height_gen.eval_index(x, y, 0));
				float const mhmin(min(min(mh00, mh01), min(mh10, mh11))), mhmax(max(max(mh00, mh01), max(mh10, mh11)));
				float const relh1(relh_adj_tex + (mhmin - zmin)*dz_inv + rand_offset);
				float const relh2(relh_adj_tex + (mhmax - zmin)*dz_inv + rand_offset);
				get_tids(relh1, NTEX_DIRT-1, h_dirt, k1, k2, t);
				get_tids(relh2, NTEX_DIRT-1, h_dirt, k3, k4, t);
				bool const same_tid(k1 == k4);
				k2 = k4;
				
				if (same_tid) {
					t = 0.0;
				}
				else {
					float const relh(relh_adj_tex + (mh00 - zmin)*dz_inv);
					get_tids(relh, NTEX_DIRT-1, h_dirt, k1, k2, t);
				}
				float const vnz(get_norm(ix).z);
				float weight_scale(1.0);
				bool const grass(lttex_dirt[k1].id == GROUND_TEX || lttex_dirt[k2].id == GROUND_TEX), snow(lttex_dirt[k2].id == SNOW_TEX);

				if (grass || snow) {
					float const *const sti(sthresh[0][snow]);

					if (vnz < sti[1]) { // handle steep slopes (dirt/rock texture replaces grass texture)
						if (grass) { // ground/grass
							float const rock_weight((lttex_dirt[k1].id == GROUND_TEX || lttex_dirt[k2].id == ROCK_TEX) ? t : 0.0);
							weight_scale = CLIP_TO_01((vnz - sti[0])/(sti[1] - sti[0]));
							weights[rock_tex_ix] += (1.0 - weight_scale)*rock_weight;
							weights[dirt_tex_ix] += (1.0 - weight_scale)*(1.0 - rock_weight);
						}
						else { // snow
							weight_scale = CLIP_TO_01(2.0f*(vnz - sti[0])/(sti[1] - sti[0]));
							weights[rock_tex_ix] += 1.0 - weight_scale;
						}
					}
				}
				weights[k2] += weight_scale*t;
				weights[k1] += weight_scale*(1.0 - t);
				unsigned const ix_val(y*tsize + x), off(4*ix_val);
				float const xv(float(x)*xy_mult), yv(float(y)*xy_mult);
				float const dirt_scale(BILINEAR_INTERP(params, dirt, xv, yv));
				float const grass_scale((mhmin < water_level) ? 0.0 : BILINEAR_INTERP(params, grass, xv, yv)); // no grass under water
				float const gscale(CLIP_TO_01(2.5f*(grass_scale - 0.5f) + 0.5f));

				if (dirt_scale  < 1.0) { // apply dirt scale: convert dirt to sand
					weights[sand_tex_ix ] += (1.0 - dirt_scale )*weights[dirt_tex_ix];
					weights[dirt_tex_ix ] *= dirt_scale;
				}
				if (grass_scale < 1.0) { // apply grass scale: convert grass to sand
					weights[sand_tex_ix ] += (1.0 - gscale)*weights[grass_tex_ix];
					weights[grass_tex_ix] *= gscale;
				}
				if (!tree_map.empty() && tree_map[ix_val] < 255) { // replace grass under trees with dirt
					float const v(tree_map[ix_val]/255.0);
					weights[dirt_tex_ix]  += (1.0 - v)*weights[grass_tex_ix];
					weights[grass_tex_ix] *= v;
				}
				for (unsigned i = 0; i < NTEX_DIRT-1; ++i) { // Note: weights should sum to 1.0, so we can calculate w4 as 1.0-w0-w1-w2-w3
					data[off+i] = (unsigned char)(255.0*CLIP_TO_01(weights[i]));
				}
				if (grass && x < size && y < size) {
					unsigned const bx(x/GRASS_BLOCK_SZ), by(y/GRASS_BLOCK_SZ), bix(by*grass_block_dim + bx);
					if (grass_blocks.empty()) {grass_blocks.resize(grass_block_dim*grass_block_dim);}
					assert(bix < grass_blocks.size());
					grass_block_t &gb(grass_blocks[bix]);
					
					if (gb.ix == 0) { // not yet set
						gb.ix   = bx + 1;
						gb.zmin = mhmin;
						gb.zmax = mhmax;
					}
					else {
						gb.zmin = min(gb.zmin, mhmin);
						gb.zmax = max(gb.zmax, mhmax);
					}
				}
			} // for x
		} // for y
		setup_texture(weight_tid, GL_MODULATE, 0, 0, 0, 0, 0);
		assert(weight_tid > 0 && glIsTexture(weight_tid));
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tsize, tsize, 0, GL_RGBA, GL_UNSIGNED_BYTE, data); // internal_format = GL_COMPRESSED_RGBA - too slow
		glDisable(GL_TEXTURE_2D);
		delete [] data;
		//PRINT_TIME("Texture Upload");
	}

	float get_rel_dist_to_camera() const {
		return max(0.0f, p2p_dist_xy(get_camera_pos(), get_center()) - radius)/get_scaled_tile_radius();
	}
	bool update_range() { // if returns 0, tile will be deleted
		if (pine_trees_enabled()) {update_pine_tree_state(0);}
		if (scenery_enabled   ()) {update_scenery();}
		if (height_tid && !grass_blocks.empty() && get_grass_dist_scale() > 1.2) {free_texture(height_tid);}
		float const dist(get_rel_dist_to_camera());
		if (dist > CLEAR_DIST_TILES) clear_vbo_tid(1,1);
		return (dist < DELETE_DIST_TILES);
	}
	bool is_visible() const {return camera_pdu.sphere_and_cube_visible_test(get_center(), radius, get_cube());}
	float get_dist_to_camera_in_tiles() const {return get_rel_dist_to_camera()*TILE_RADIUS;}
	float get_scenery_dist_scale () const {return get_dist_to_camera_in_tiles()/SCENERY_THRESH;}
	float get_tree_dist_scale    () const {return get_dist_to_camera_in_tiles()/max(1.0f, TREE_LOD_THRESH*calc_tree_size());}
	float get_tree_far_weight    () const {return (ENABLE_TREE_LOD ? CLIP_TO_01(GEOMORPH_THRESH*(get_tree_dist_scale() - 1.0f)) : 0.0);}
	float get_grass_dist_scale   () const {return get_dist_to_camera_in_tiles()/GRASS_THRESH;}

	// *** pine trees ***

	template <typename T> void postproc_trees(T const &trees, float &tzmax) { // pine/decidious trees
		tzmax  = mzmin;
		trees.update_zmax(tzmax);
		trmax  = max(trmax, trees.get_rmax());
		radius = max(radius, (calc_radius() + trmax)); // is this really needed?
	}

	void init_pine_tree_draw() {
		ptree_off.set_from_xyoff2();
		pine_trees.gen_trees(x1+ptree_off.dxoff, y1+ptree_off.dyoff, x2+ptree_off.dxoff, y2+ptree_off.dyoff, vegetation*get_avg_veg());
		pine_trees.calc_trunk_pts();
		postproc_trees(pine_trees, ptzmax);
	}

	void update_pine_tree_state(bool upload_if_needed) {
		if (pine_trees.empty()) return;
		float const weight(get_tree_far_weight());
		float const weights[2] = {1.0-weight, weight}; // {high, low} detail

		for (unsigned d = 0; d < 2; ++d) {
			if (weights[d] > 0.0) { // needed
				if (upload_if_needed) {
					vector3d const delta(get_camera_pos() - get_center());
					bool const pri_dim(fabs(delta.y) < fabs(delta.x));
					pine_trees.finalize_upload_and_clear_pts((d != 0), pri_dim); // needed for drawing
				}
			}
			else { // not needed
				pine_trees.clear_vbo_and_ids_if_needed(d != 0);
			}
		}
	}

	unsigned num_pine_trees() const {return pine_trees.size();}

	void draw_tree_leaves_lod(shader_t &s, vector3d const &xlate, bool low_detail) const {
		bool const draw_all(low_detail || camera_pdu.point_visible_test(get_center())); // tile center is in view
		pine_trees.draw_pine_leaves(0, low_detail, draw_all, xlate);
	}

	void draw_pine_trees(shader_t &s, vector<point> &trunk_pts, bool draw_branches, bool draw_near_leaves, bool draw_far_leaves, bool reflection_pass) const {
		if (pine_trees.empty()) return;
		glPushMatrix();
		vector3d const xlate(ptree_off.get_xlate());
		translate_to(xlate);
		
		if (draw_branches) {
			float const dscale(get_tree_dist_scale());

			if (dscale < 1.0) { // close, draw as polygons
				pine_trees.draw_branches(0, xlate, &trunk_pts);
			}
			else if (dscale < 2.0) { // far away, use low detail branches
				pine_trees.add_trunk_pts(xlate, trunk_pts);
			} // else very far, skip branches
		}
		if (draw_near_leaves || draw_far_leaves) { // could use reflection_pass as an optimization
			float const weight(1.0 - get_tree_far_weight()); // 0 => low detail, 1 => high detail

			if (weight > 0 && weight < 1.0) { // use geomorphing with dithering (since alpha doesn't blend in the correct order)
				if (draw_near_leaves) { // GL_SAMPLE_ALPHA_TO_COVERAGE?
					s.add_uniform_float("max_noise", weight);
					draw_tree_leaves_lod(s, xlate, 0);
					s.add_uniform_float("max_noise", 1.0);
				}
				if (draw_far_leaves) {
					s.add_uniform_float("min_noise", weight);
					draw_tree_leaves_lod(s, xlate, 1);
					s.add_uniform_float("min_noise", 0.0);
				}
			}
			else if ((weight == 0.0) ? draw_far_leaves : draw_near_leaves) {
				draw_tree_leaves_lod(s, xlate, (weight == 0.0));
			}
		}
		glPopMatrix();
	}

	// *** decidious trees ***

	unsigned num_decid_trees() const {return decid_trees.size();}

	void gen_decid_trees_if_needed() {
		if (decid_trees.was_generated()) return; // already generated
		assert(decid_trees.empty());
		dtree_off.set_from_xyoff2();
		decid_trees.gen_deterministic(x1+dtree_off.dxoff, y1+dtree_off.dyoff, x2+dtree_off.dxoff, y2+dtree_off.dyoff, vegetation*get_avg_veg());
		postproc_trees(decid_trees, dtzmax);
	}

	void draw_decid_trees(shader_t &s, bool draw_branches, bool draw_leaves, bool reflection_pass) {
		if (decid_trees.empty()) return;
		decid_trees.draw_branches_and_leaves(s, draw_branches, draw_leaves, 0, dtree_off.get_xlate());
	}

	// *** scenery ***

	void update_scenery() {
		float const dist_scale(get_scenery_dist_scale()); // tree_dist_scale should correlate with mesh scale
		if (scenery.generated && dist_scale > 1.2) {scenery.clear();} // too far away
		if (scenery.generated || dist_scale > 1.0 || !is_visible()) return; // already generated, too far away, or not visible
		scenery_off.set_from_xyoff2();
		scenery.gen(x1+scenery_off.dxoff, y1+scenery_off.dyoff, x2+scenery_off.dxoff, y2+scenery_off.dyoff, vegetation*get_avg_veg());
	}

	void draw_scenery(shader_t &s, bool draw_opaque, bool draw_leaves) {
		if (!scenery.generated || get_scenery_dist_scale() > 1.0) return;
		glPushMatrix();
		vector3d const xlate(scenery_off.get_xlate());
		translate_to(xlate);
		if (draw_opaque) {scenery.draw_opaque_objects(s, 0, xlate, 1);}
		if (draw_leaves) {scenery.draw_plant_leaves  (s, 0, xlate);}
		glPopMatrix();
	}

	// *** grass ***

	void init_draw_grass() {
		if (height_tid || grass_blocks.empty() || get_grass_dist_scale() > 1.0) return;
		float const scale(65535/(mzmax - mzmin));
		assert(zvals.size() == zvsize*zvsize);
		vector<unsigned short> data(stride*stride);

		for (unsigned y = 0; y < stride; ++y) {
			for (unsigned x = 0; x < stride; ++x) {
				data[y*stride+x] = (unsigned short)(scale*(zvals[y*zvsize+x] - mzmin));
			}
		}
		setup_texture(height_tid, GL_MODULATE, 0, 0, 0, 0, 0);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 2);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE16, stride, stride, 0, GL_LUMINANCE, GL_UNSIGNED_SHORT, &data.front());
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	}

	void draw_grass(shader_t &s, bool use_cloud_shadows) {
		if (grass_blocks.empty() || get_grass_dist_scale() > 1.0) return;
		bind_texture_tu(height_tid, 2);
		bind_texture_tu(weight_tid, 3);
		bind_texture_tu(shadow_normal_tid, 4);
		s.add_uniform_float("x1",  -0.5*DX_VAL);
		s.add_uniform_float("y1",  -0.5*DY_VAL);
		s.add_uniform_float("x2",  (size + 0.5)*DX_VAL);
		s.add_uniform_float("y2",  (size + 0.5)*DY_VAL);
		s.add_uniform_float("wx2", size*DX_VAL);
		s.add_uniform_float("wy2", size*DY_VAL);
		s.add_uniform_float("zmin", mzmin);
		s.add_uniform_float("zmax", mzmax);
		
		if (use_cloud_shadows) {
			vector3d const offset(get_xval(x1), get_yval(y1), 0.0);
			s.add_uniform_vector3d("cloud_offset", offset);
		}
		unsigned const grass_block_dim(get_grass_block_dim());
		assert(grass_blocks.size() == grass_block_dim*grass_block_dim);
		unsigned const ty_loc(s.get_uniform_loc("translate_y"));
		int const dx(xoff - xoff2), dy(yoff - yoff2);
		float const llcx(get_xval(x1+dx)), llcy(get_yval(y1+dy)), dx_step(GRASS_BLOCK_SZ*DX_VAL), dy_step(GRASS_BLOCK_SZ*DY_VAL);
		float const lod_scale(1.0/get_scaled_tile_radius());
		point const camera(get_camera_pos());
		glPushMatrix();
		glTranslatef(llcx, llcy, 0.0);

		for (unsigned y = 0; y < grass_block_dim; ++y) {
			s.set_uniform_float(ty_loc, y*dy_step);

			for (unsigned x = 0; x < grass_block_dim; ++x) {
				grass_block_t const &gb(grass_blocks[y*grass_block_dim+x]);
				if (gb.ix == 0) continue; // empty block
				cube_t const bcube(llcx+x*dx_step, llcx+(x+1)*dx_step, llcy+y*dy_step, llcy+(y+1)*dy_step, gb.zmin, (gb.zmax + grass_length));
				point const center(bcube.get_cube_center());
				if (max(0.0f, p2p_dist_xy(camera, center) - radius)*TILE_RADIUS*lod_scale > GRASS_THRESH) continue;
				if (!camera_pdu.cube_visible(bcube)) continue;
				bool back_facing(1);

				for (unsigned yy = y*GRASS_BLOCK_SZ; yy <= (y+1)*GRASS_BLOCK_SZ && back_facing; ++yy) {
					for (unsigned xx = x*GRASS_BLOCK_SZ; xx <= (x+1)*GRASS_BLOCK_SZ && back_facing; ++xx) {
						unsigned const ix(yy*zvsize + xx);
						back_facing &= (dot_product(get_norm(ix), (camera - point(llcx+xx*DX_VAL, llcy+yy*DY_VAL, zvals[ix]))) < 0.0);
					}
				}
				if (back_facing) continue;
				unsigned const lod_level(min(NUM_GRASS_LODS-1, unsigned(GRASS_LOD_SCALE*lod_scale*distance_to_camera(center))));
				grass_tile_manager.render_block((gb.ix - 1), lod_level);
			}
		}
		glPopMatrix();
	}

	// *** rendering ***

	void ensure_vbo(vector<vert_type_t> &data, vector<unsigned short> indices[NUM_LODS]) {
		if (vbo) return; // already allocated
		create_data(data, indices);
		if (any_trees_enabled()) {apply_tree_ao_shadows();}
		create_vbo_and_upload(vbo, data, 0, 0);

		for (unsigned i = 0; i < NUM_LODS; ++i) {
			assert(ivbo[i] == 0);
			assert(!indices[i].empty());
			create_vbo_and_upload(ivbo[i], indices[i], 1, 0);
		}
	}

	void ensure_weights(mesh_xy_grid_cache_t &height_gen) {
		if (weights_invalid) {free_texture(weight_tid);}
		if (weight_tid == 0) {create_texture(height_gen);}
		check_shadow_map_and_normal_texture();
	}

	void draw(shader_t &s, bool reflection_pass) const {
		assert(vbo);
		assert(size > 0);
		assert(weight_tid > 0);
		bind_2d_texture(weight_tid);
		glPushMatrix();
		translate_to(mesh_off.get_xlate());
		set_landscape_texgen(1.0, (-x1 - mesh_off.dxoff), (-y1 - mesh_off.dyoff), MESH_X_SIZE, MESH_Y_SIZE);
		
		if (!reflection_pass && cloud_shadows_enabled()) {
			vector3d const offset(-mesh_off.dxoff*DX_VAL, -mesh_off.dyoff*DY_VAL, 0.0);
			s.add_uniform_vector3d("cloud_offset", offset);
		}
		bind_texture_tu(shadow_normal_tid, 7);
		unsigned lod_level(reflection_pass ? min(NUM_LODS-1, 1U) : 0);
		float dist(get_dist_to_camera_in_tiles());

        while (dist > (reflection_pass ? 1.0 : 2.0) && lod_level+1 < NUM_LODS) {
            dist /= 2;
            ++lod_level;
        }
		assert(lod_level < NUM_LODS);
		assert(vbo > 0 && ivbo[lod_level] > 0);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo[lod_level], 1);
		unsigned const isz(size >> lod_level), ptr_stride(sizeof(vert_type_t));
		
		// can store normals in a normal map texture, but a vertex texture fetch is slow
		glVertexPointer(3, GL_FLOAT, ptr_stride, 0);
		glDrawRangeElements(GL_QUADS, 0, stride*stride, 4*isz*isz, GL_UNSIGNED_SHORT, 0); // requires GL/glew.h
		bind_vbo(0, 0);
		bind_vbo(0, 1);
		glPopMatrix();
		if (weight_tid > 0) disable_textures_texgen();
	}

	void draw_water(float z) const {
		if (!has_water() || get_rel_dist_to_camera() > DRAW_DIST_TILES || !is_visible()) return;
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2)), xv2(xv1+(x2-x1)*DX_VAL), yv2(yv1+(y2-y1)*DY_VAL);
		glVertex3f(xv1, yv1, z);
		glVertex3f(xv1, yv2, z);
		glVertex3f(xv2, yv2, z);
		glVertex3f(xv2, yv1, z);
	}
}; // tile_t


class tile_draw_t {

	typedef map<tile_xy_pair, tile_t*> tile_map;
	tile_map tiles;
	vector<point> tree_trunk_pts;
	mesh_xy_grid_cache_t height_gen;

public:
	tile_draw_t() {
		assert(MESH_X_SIZE == MESH_Y_SIZE && X_SCENE_SIZE == Y_SCENE_SIZE);
	}
	//~tile_draw_t() {clear();}

	void clear() {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->clear();
			delete i->second;
		}
		tiles.clear();
		tree_trunk_pts.clear();
	}

	void update() {
		//RESET_TIME;
		grass_tile_manager.update(); // every frame, even if not in tiled terrain mode?
		assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
		point const camera(get_camera_pos() - point((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
		int const tile_radius(int(1.5*TILE_RADIUS) + 1);
		int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
		int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
		int const x2( tile_radius + toffx), y2( tile_radius + toffy);
		unsigned const init_tiles((unsigned)tiles.size());
		vector<tile_xy_pair> to_erase;

		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) { // update tiles and free old tiles
			if (!i->second->update_range()) {
				to_erase.push_back(i->first);
				i->second->clear();
				delete i->second;
			}
		}
		for (vector<tile_xy_pair>::const_iterator i = to_erase.begin(); i != to_erase.end(); ++i) {
			tiles.erase(*i);
		}
		for (int y = y1; y <= y2; ++y ) { // create new tiles
			for (int x = x1; x <= x2; ++x ) {
				tile_xy_pair const txy(x, y);
				if (tiles.find(txy) != tiles.end()) continue; // already exists
				tile_t tile(MESH_X_SIZE, x, y);
				if (tile.get_rel_dist_to_camera() >= CREATE_DIST_TILES) continue; // too far away to create
				tile_t *new_tile(new tile_t(tile));
				new_tile->create_zvals(height_gen);
				tiles[txy] = new_tile;
			}
		}
		if (DEBUG_TILES && (tiles.size() != init_tiles || !to_erase.empty())) {
			cout << "update: tiles: " << init_tiles << " to " << tiles.size() << ", erased: " << to_erase.size() << endl;
		}
		//PRINT_TIME("Tiled Terrain Update");
	}


	static void setup_terrain_textures(shader_t &s, unsigned start_tu_id, bool use_sand) {
		unsigned const base_tsize(NORM_TEXELS);

		for (int i = 0; i < (use_sand ? NTEX_SAND : NTEX_DIRT); ++i) {
			int const tid(use_sand ? lttex_sand[i].id : lttex_dirt[i].id);
			float const tscale(float(base_tsize)/float(get_texture_size(tid, 0))); // assumes textures are square
			float const cscale((tid == GROUND_TEX) ? TT_GRASS_COLOR_SCALE : 1.0); // darker grass
			unsigned const tu_id(start_tu_id + i);
			select_multitex(tid, tu_id, 0);
			std::ostringstream oss1, oss2, oss3;
			oss1 << "tex" << tu_id;
			oss2 << "ts"  << tu_id;
			oss3 << "cs"  << tu_id;
			s.add_uniform_int(  oss1.str().c_str(), tu_id);
			s.add_uniform_float(oss2.str().c_str(), tscale);
			s.add_uniform_float(oss3.str().c_str(), cscale);
		}
		set_multitex(0);
	}

	// uses texture units 0-9
	static void setup_mesh_draw_shaders(shader_t &s, bool reflection_pass) {
		s.setup_enabled_lights();
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		s.set_prefix("#define NUM_OCTAVES 4",     1); // FS (for clouds)
		s.set_bool_prefix("apply_cloud_shadows", (cloud_shadows_enabled() && !reflection_pass), 1); // FS
		s.set_vert_shader("texture_gen.part+water_fog.part+tiled_mesh");
		s.set_frag_shader("linear_fog.part+perlin_clouds.part*+tiled_mesh");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		s.add_uniform_int("tex1", 1);
		s.add_uniform_int("shadow_normal_tex", 7);
		s.add_uniform_float("water_plane_z", get_tiled_terrain_water_level());
		s.add_uniform_float("water_atten", WATER_COL_ATTEN*mesh_scale);
		s.add_uniform_float("normal_z_scale", (reflection_pass ? -1.0 : 1.0));
		set_noise_tex(s, 8);
		s.add_uniform_float("cloud_scale", 0.535);
		s.add_uniform_float("cloud_plane_z", get_cloud_zmax());
		set_cloud_uniforms(s, 9);
		setup_terrain_textures(s, 2, 0);
	}

	typedef vector<pair<float, tile_t *> > draw_vect_t;

	float draw(bool reflection_pass) {
		float zmin(FAR_CLIP);
		set_array_client_state(1, 0, 0, 0);
		unsigned num_drawn(0), num_trees(0);
		unsigned long long mem(grass_tile_manager.get_gpu_mem() + tree_data_manager.get_gpu_mem()), tree_mem(0);
		vector<tile_t::vert_type_t> data;
		vector<unsigned short> indices[NUM_LODS];
		static point last_sun(all_zeros), last_moon(all_zeros);
		draw_vect_t to_draw;
		vector<tile_t *> to_gen_trees;

		if (sun_pos != last_sun || moon_pos != last_moon) { // light source change
			for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
				i->second->clear_shadows();
			}
			last_sun  = sun_pos;
			last_moon = moon_pos;
		}
		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
			tile_t *const tile(i->second);
			assert(tile);
			
			if (DEBUG_TILES) {
				mem      += tile->get_gpu_mem();
				tree_mem += tile->get_tree_mem();
			}
			float const dist(tile->get_rel_dist_to_camera());
			if (dist > DRAW_DIST_TILES) continue; // too far to draw
			zmin = min(zmin, tile->get_zmin());
			if (!tile->is_visible())    continue;
			if (pine_trees_enabled() && !tile->pine_trees_generated()) {to_gen_trees.push_back(tile);}
			if (decid_trees_enabled()) {tile->gen_decid_trees_if_needed();}
			if (reflection_pass && ((tile->contains_camera() && !tile->has_water()) || tile->all_water())) continue;
			to_draw.push_back(make_pair(dist, tile));
		}
		if (!to_gen_trees.empty()) {
			#pragma omp parallel for schedule(static,1)
			for (int i = 0; i < (int)to_gen_trees.size(); ++i) {
				to_gen_trees[i]->init_pine_tree_draw();
			}
		}
		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->ensure_vbo(data, indices);
		}
		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->ensure_weights(height_gen);
		}
		sort(to_draw.begin(), to_draw.end()); // sort front to back to improve draw time through depth culling
		shader_t s;
		setup_mesh_draw_shaders(s, reflection_pass);
		s.add_uniform_float("spec_scale", 1.0);
		setup_mesh_lighting();

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			num_trees += to_draw[i].second->num_pine_trees();
			if (display_mode & 0x01) {to_draw[i].second->draw(s, reflection_pass);}
		}
		s.end_shader();
		
		if (DEBUG_TILES) {
			cout << "tiles drawn: " << to_draw.size() << " of " << tiles.size() << ", trees drawn: "
				<< num_trees << ", gpu mem: " << mem/1024/1024 << ", tree mem: " << tree_mem/1024/1024 << endl;
		}
		run_post_mesh_draw();
		if (pine_trees_enabled ()) {draw_pine_trees (to_draw, reflection_pass);}
		if (decid_trees_enabled()) {draw_decid_trees(to_draw, reflection_pass);}
		if (scenery_enabled    ()) {draw_scenery    (to_draw, reflection_pass);}
		if (is_grass_enabled   ()) {draw_grass      (to_draw, reflection_pass);}
		return zmin;
	}

	void draw_water(float zval) const {
		glBegin(GL_QUADS);

		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->draw_water(zval);
		}
		glEnd();
	}

	static void set_noise_tex(shader_t &s, unsigned tu_id) {
		set_multitex(tu_id);
		select_texture(NOISE_GEN_TEX, 0);
		s.add_uniform_int("noise_tex", tu_id);
		set_multitex(0);
	}

	void set_pine_tree_shader(shader_t &s, string const &vs) const {
		s.set_prefix("#define USE_LIGHT_COLORS",   0); // VS
		s.set_prefix("#define USE_GOOD_SPECULAR",  0); // VS
		s.set_prefix("#define USE_QUADRATIC_FOG",  1); // FS
		s.setup_enabled_lights(2);
		s.set_vert_shader("ads_lighting.part*+tc_by_vert_id.part+" + vs);
		s.set_frag_shader("linear_fog.part+pine_tree");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("branch_tex", 0);
		s.add_uniform_float("min_alpha", 0.75);

		set_noise_tex(s, 1);
		s.add_uniform_float("noise_tex_size", get_texture_size(NOISE_GEN_TEX, 0));
		check_gl_error(302);
	}

	void draw_pine_tree_bl(draw_vect_t const &to_draw, shader_t &s, bool branches, bool near_leaves, bool far_leaves, bool reflection_pass) {
		for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
			to_draw[i].second->draw_pine_trees(s, tree_trunk_pts, branches, near_leaves, far_leaves, reflection_pass);
		}
	}

	void draw_pine_trees(draw_vect_t const &to_draw, bool reflection_pass) {
		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->update_pine_tree_state(1);
		}

		// nearby trunks
		shader_t s;
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0);
		s.add_uniform_color("const_indir_color", colorRGB(0,0,0)); // don't want indir lighting for tree trunks
		s.add_uniform_float("tex_scale_t", 5.0);
		set_color(get_tree_trunk_color(T_PINE, 0)); // all a constant color
		draw_pine_tree_bl(to_draw, s, 1, 0, 0, reflection_pass); // branches
		s.add_uniform_float("tex_scale_t", 1.0);
		s.end_shader();

		// leaves/distant trunks
		set_pine_tree_shader(s, "pine_tree_billboard");
		
		if (!tree_trunk_pts.empty()) { // color/texture already set above
			assert(!(tree_trunk_pts.size() & 1));
			select_texture(WHITE_TEX, 0); // enable=0
			get_tree_trunk_color(T_PINE, 1).do_glColor();
			zero_vector.do_glNormal();
			set_array_client_state(1, 0, 0, 0);
			glVertexPointer(3, GL_FLOAT, sizeof(point), &tree_trunk_pts.front());
			glDrawArrays(GL_LINES, 0, (unsigned)tree_trunk_pts.size());
			tree_trunk_pts.resize(0);
		}
		set_specular(0.2, 8.0);
		draw_pine_tree_bl(to_draw, s, 0, 0, 1, reflection_pass); // far leaves
		s.end_shader();
		set_pine_tree_shader(s, "pine_tree");
		draw_pine_tree_bl(to_draw, s, 0, 1, 0, reflection_pass); // near leaves
		assert(tree_trunk_pts.empty());
		s.end_shader();
		set_specular(0.0, 1.0);
	}

	void draw_decid_tree_bl(draw_vect_t const &to_draw, shader_t &s, bool branches, bool leaves, bool reflection_pass) {
		for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
			to_draw[i].second->draw_decid_trees(s, branches, leaves, reflection_pass);
		}
	}

	void draw_decid_trees(draw_vect_t const &to_draw, bool reflection_pass) {
		//unsigned tot(0);
		//for (unsigned i = 0; i < to_draw.size(); ++i) {tot += to_draw[i].second->num_decid_trees();}
		//cout << "to draw: " << to_draw.size() << " of " << tiles.size() << ", total trees: " << tot << endl;

		// FIXME: faster (view clipping, LOD, etc.)
		// FIXME: shadow map?

		// draw branches
		shader_t bs;
		bs.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		set_tree_branch_shader(bs, 1, 0, 0);
		bs.add_uniform_color("const_indir_color", colorRGB(0,0,0)); // don't want indir lighting for tree trunks
		draw_decid_tree_bl(to_draw, bs, 1, 0, reflection_pass);
		bs.add_uniform_vector3d("world_space_offset", zero_vector); // reset
		bs.end_shader();
		disable_multitex_a();

		// draw leaves
		shader_t ls;
		tree_cont_t::pre_leaf_draw(ls);
		float const cscale(cloud_shadows_enabled() ? 0.75 : 1.0);
		ls.add_uniform_color("color_scale", colorRGBA(cscale, cscale, cscale, 1.0));
		draw_decid_tree_bl(to_draw, ls, 0, 1, reflection_pass);
		ls.add_uniform_color("color_scale", WHITE);
		tree_cont_t::post_leaf_draw(ls);
		leaf_color_changed = 0; // Note: only visible trees will be updated
	}

	void draw_scenery(draw_vect_t const &to_draw, bool reflection_pass) {
		shader_t s;
		s.setup_enabled_lights(2);
		s.set_prefix("#define USE_LIGHT_COLORS",  0); // VS
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		s.set_vert_shader("ads_lighting.part*+two_lights_texture");
		s.set_frag_shader("linear_fog.part+textured_with_fog");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		
		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->draw_scenery(s, 1, 0); // opaque
		}
		s.end_shader();
		set_leaf_shader(s, 0.9, 1, 0);

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->draw_scenery(s, 0, 1); // leaves
		}
		s.end_shader();
	}

	// tu's used: 0: grass, 1: wind noise, 2: heightmap, 3: grass weight, 4: shadow map, 5: noise, 9: cloud noise
	void draw_grass(draw_vect_t const &to_draw, bool reflection_pass) {
		if (reflection_pass) return; // no grass refletion (yet)

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			to_draw[i].second->init_draw_grass();
		}
		grass_tile_manager.begin_draw(0.1);
		shader_t s;
		bool const use_cloud_shadows(GRASS_CLOUD_SHADOWS && cloud_shadows_enabled() && !reflection_pass);

		for (unsigned pass = 0; pass < 2; ++pass) { // wind, no wind
			s.set_prefix("#define USE_LIGHT_COLORS",  0); // VS
			s.set_prefix("#define USE_GOOD_SPECULAR", 0); // VS
			s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
			s.set_prefix("#define NUM_OCTAVES 4",     0); // VS (for clouds)
			s.set_bool_prefix("apply_cloud_shadows", use_cloud_shadows, 0); // VS
			s.set_bool_prefix("enable_grass_wind", (pass == 0), 0); // VS
			s.setup_enabled_lights(2);
			s.set_vert_shader("ads_lighting.part*+wind.part*+perlin_clouds.part*+grass_tiled");
			s.set_frag_shader("linear_fog.part+textured_with_fog");
			//s.set_geom_shader("ads_lighting.part*+grass_tiled", GL_TRIANGLES, GL_TRIANGLE_STRIP, 3); // too slow
			s.begin_shader();
			if (pass == 0) {setup_wind_for_shader(s, 1);}
			s.add_uniform_int("tex0", 0);
			s.add_uniform_int("height_tex", 2);
			s.add_uniform_int("weight_tex", 3);
			s.add_uniform_int("shadow_normal_tex", 4);
			set_noise_tex(s, 5);
			s.setup_fog_scale();
			s.add_uniform_float("height", grass_length);
			s.add_uniform_float("dist_const", (X_SCENE_SIZE + Y_SCENE_SIZE)*GRASS_THRESH);
			s.add_uniform_float("dist_slope", 0.25);
			s.add_uniform_float("cloud_scale", 0.535);
			s.add_uniform_float("cloud_plane_z", get_cloud_zmax());
			set_cloud_uniforms(s, 9);

			for (unsigned i = 0; i < to_draw.size(); ++i) {
				if ((to_draw[i].second->get_dist_to_camera_in_tiles() > 0.5) == pass) {
					to_draw[i].second->draw_grass(s, use_cloud_shadows);
				}
			}
			s.end_shader();
		}
		grass_tile_manager.end_draw();
	}

	void clear_vbos_tids(bool vclear, bool tclear) {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->clear_vbo_tid(vclear, tclear);
		}
	}

	tile_t *get_tile_from_xy(tile_xy_pair const &tp) {
		tile_map::iterator it(tiles.find(tp));
		if (it != tiles.end()) return it->second;
		return NULL;
	}
}; // tile_draw_t


tile_draw_t terrain_tile_draw;


tile_t *get_tile_from_xy(tile_xy_pair const &tp) {

	return terrain_tile_draw.get_tile_from_xy(tp);
}


float draw_tiled_terrain(bool reflection_pass) {

	//RESET_TIME;
	bool const vbo_supported(setup_gen_buffers());
		
	if (!vbo_supported) {
		cout << "Warning: VBOs not supported, so tiled mesh cannot be drawn." << endl;
		return zmin;
	}
	terrain_tile_draw.update();
	float const zmin(terrain_tile_draw.draw(reflection_pass));
	//glFinish(); PRINT_TIME("Tiled Terrain Draw"); //exit(0);
	return zmin;
}


void clear_tiled_terrain() {
	terrain_tile_draw.clear();
}

void reset_tiled_terrain_state() {
	terrain_tile_draw.clear_vbos_tids(1,1);
}

void draw_tiled_terrain_water(float zval) {
	terrain_tile_draw.draw_water(zval);
}


