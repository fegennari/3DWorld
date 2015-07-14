// 3D World - Tiled Landscape Mesh Generation/Drawing Code
// by Frank Gennari
// 9/26/10

#include "tiled_mesh.h"
#include "textures_3dw.h"
#include "shaders.h"
#include "openal_wrap.h"
#include "heightmap.h"


bool const DEBUG_TILES        = 0;
bool const DEBUG_TILE_BOUNDS  = 0;
bool const ENABLE_INST_PINE   = 1; // faster generation, lower GPU memory, slower rendering
int  const DITHER_NOISE_TEX   = NOISE_GEN_TEX;//PS_NOISE_TEX
unsigned const NORM_TEXELS    = 512;
unsigned const NUM_FIRE_MODES = 4;
float const FOG_DIST_TILES    = 1.45;
float const DRAW_DIST_TILES   = 1.5;
float const CREATE_DIST_TILES = 1.6;
float const CLEAR_DIST_TILES  = 1.6;
float const DELETE_DIST_TILES = 1.8;
float const GRASS_LOD_SCALE   = 15.0; // smaller = more grass detail
float const GRASS_DIST_SLOPE  = 0.25;
float const SMAP_FADE_THRESH  = 1.5;
float const OCCLUDER_DIST     = 0.2;
float const FLOWER_REL_DIST   = 0.9; // flower view distance relative to grass view distance

int   const LIGHTNING_LIGHT = 2;
float const LIGHTNING_FREQ  = 200.0; // in ticks (1/40 s)
float const LITNING_TIME2   = 40.0;
float const LITNING_DIST    = 1.2;


bool tt_lightning_enabled(0);
unsigned inf_terrain_fire_mode(0); // none, increase height, decrease height
string read_hmap_modmap_fn, write_hmap_modmap_fn("heightmap.mod");
hmap_brush_param_t cur_brush_param;
tile_t::offset_t model3d_offset;

extern bool inf_terrain_scenery, enable_tiled_mesh_ao, underwater, fog_enabled, volume_lighting;
extern unsigned grass_density, max_unique_trees, inf_terrain_fire_mode, shadow_map_sz;
extern int DISABLE_WATER, display_mode, tree_mode, leaf_color_changed, ground_effects_level, animate2, iticks, num_trees;
extern int invert_mh_image, is_cloudy, camera_surf_collide, show_fog, mesh_gen_mode, mesh_gen_shape, cloud_model;
extern float zmax, zmin, water_plane_z, mesh_scale, mesh_scale_z, vegetation, relh_adj_tex, grass_length, grass_width, fticks, tfticks;
extern float ocean_wave_height, sm_tree_density, tree_scale, atmosphere, cloud_cover, temperature, flower_density, FAR_CLIP, shadow_map_pcf_offset;
extern point sun_pos, moon_pos, surface_pos;
extern vector3d wind;
extern water_params_t water_params;
extern char *mh_filename_tt;
extern float h_dirt[];
extern tree_data_manager_t tree_data_manager;
extern pt_line_drawer tree_scenery_pld;

bool enable_terrain_env(ENABLE_TERRAIN_ENV);
void set_water_plane_uniforms(shader_t &s);

void create_pine_tree_instances();
unsigned get_pine_tree_inst_gpu_mem();

void setup_detail_normal_map(shader_t &s, float tscale);
void draw_distant_mesh_bottom(float terrain_zmin);


float get_inf_terrain_fog_dist() {return FOG_DIST_TILES*get_scaled_tile_radius();}
float get_draw_tile_dist  () {return DRAW_DIST_TILES*get_scaled_tile_radius();}
float get_grass_thresh    () {return GRASS_THRESH*get_tile_width();}
float get_grass_thresh_pad() {return (get_grass_thresh() + 1.0/GRASS_DIST_SLOPE);}
bool is_water_enabled     () {return (!DISABLE_WATER && (display_mode & 0x04) != 0);}
bool pine_trees_enabled   () {return ((tree_mode & 2) && vegetation > 0.0);}
bool decid_trees_enabled  () {return ((tree_mode & 1) && vegetation > 0.0);}
bool any_trees_enabled    () {return (pine_trees_enabled() || decid_trees_enabled());}
bool scenery_enabled      () {return (inf_terrain_scenery && SCENERY_THRESH > 0.0);}
bool gen_grass_map        () {return (GRASS_THRESH > 0.0 && grass_density > 0 && vegetation > 0.0);}
bool is_grass_enabled     () {return ((display_mode & 0x02) && gen_grass_map());}
bool clouds_enabled       () {return ((display_mode & 0x40) == 0);}
bool cloud_shadows_enabled() {return (ground_effects_level >= 2 && clouds_enabled());}
bool mesh_shadows_enabled () {return (ground_effects_level >= 1);}
bool is_distance_mode     () {return ((display_mode & 0x10) != 0);}
bool nonunif_fog_enabled  () {return (show_fog && is_distance_mode());}
bool enable_ocean_waves   () {return ((display_mode & 0x0100) != 0 && wind.mag() > TOLERANCE);}
bool draw_distant_water   () {return (is_water_enabled() && is_distance_mode() && camera_pdu.far_ > 1.1*FAR_CLIP);}
bool use_water_plane_tess () {return (enable_ocean_waves() && cloud_model == 0 && !draw_distant_water());} // hack to use cloud_model (F10)
float get_tt_fog_top      () {return (nonunif_fog_enabled() ? (zmax + (zmax - zmin)) : (zmax + FAR_CLIP));}
float get_tt_fog_bot      () {return (nonunif_fog_enabled() ? zmax : (zmax + FAR_CLIP));}
float get_tt_cloud_level  () {return 0.5*(get_tt_fog_bot() + get_tt_fog_top());}
float get_smap_atten_val  () {return SMAP_FADE_THRESH*get_tile_width();}
unsigned get_tile_size    () {return MESH_X_SIZE;}

bool enable_instanced_pine_trees() {
	float const ntrees_mult(vegetation*sm_tree_density*tree_scale*tree_scale);
	return (ENABLE_INST_PINE && ntrees_mult > 20 && max_unique_trees > 0); // enable when there are lots of trees
}


float get_tt_fog_based_far_clip(float min_camera_dist) {

	float const uniform_fog_far_clip(FAR_CLIP + 2.0*min_camera_dist);
	if (!nonunif_fog_enabled()) {return uniform_fog_far_clip;}
	float const zf(get_tt_fog_bot()), z0(get_tt_fog_top()), zc(get_camera_pos().z), zm(zmax);
	assert(zm <= zf && zf <= z0); // we currently only support these cases
	if (zc < zf) {return uniform_fog_far_clip;} // camera below max fog line, fog density == 1.0, so use standard far clip
	float const FD(get_inf_terrain_fog_dist()); // 67.2
	float dfavg;
	
	if (zc < z0) { // camera within linear fog/atmosphere region
		dfavg = (2*z0 - zf - zc)/(2*(z0 - zf))*(zc - zf) + (zf - zm);
	}
	else { // camera above fog/atmosphere
		dfavg = 0.5*(zf + z0) - zm;
	}
	float const fog_dist((zc - zm)*sqrt(max((FD*FD/(dfavg*dfavg) - 1.0), 0.0)));
	//cout << "FD: " << FD << ", fog_dist: " << fog_dist << ", FAR_CLIP: " << FAR_CLIP << ", final: " << max(FAR_CLIP, 1.1f*fog_dist) << endl;
	return max(uniform_fog_far_clip, 1.1f*fog_dist); // add 10% padding
}


grass_tile_manager_t grass_tile_manager;

void update_tiled_terrain_grass_vbos() {grass_tile_manager.clear_vbo();}
void update_tiled_grass_colors() {grass_tile_manager.clear();} // regenerate grass


void bind_texture_tu(unsigned tid, unsigned tu_id) {

	assert(tid);
	set_active_texture(tu_id);
	bind_2d_texture(tid);
	set_active_texture(0);
}


#define BILINEAR_INTERP(arr, var, x, y) (y*(x*arr[1][1].var + (1.0-x)*arr[1][0].var) + (1.0-y)*(x*arr[0][1].var + (1.0-x)*arr[0][0].var))


// *** heightmap management ***


class tiled_terrain_hmap_manager_t : public terrain_hmap_manager_t {

	tile_t *cur_tile;
	bool modified[3][3];

public:
	tiled_terrain_hmap_manager_t() : cur_tile(NULL) {clear_modified();}
	void clear_modified() {for (unsigned i = 0; i < 3; ++i) {UNROLL_3X(modified[i][i_] = 0;)}}

	void apply_brush(tex_mod_map_manager_t::hmap_brush_t brush, tile_t *tile, bool cache) { // Note: brush is copied and may be modified
		cur_tile = tile;
		assert(brush.radius <= get_tile_size()); // only allow for a single adjacent tile
		clear_modified();

		if (brush.is_flatten_brush()) { // use heightmap value at brush center instead of a delta
			brush.delta = get_clamped_pixel_value(brush.x, brush.y); // Note: original delta is overwritten/unused in this case
		}
		int const step_sz(max(1, int(1.0/mesh_scale + SMALL_NUMBER))); // Note: only intended to work when mesh_scale is a power of 0.5 (or generally an integer reciprocol)
		unsigned const num_steps(max(1U, unsigned(mesh_scale + SMALL_NUMBER))); // Note: only intended to work when mesh_scale is a power of 2 (or generally an integer)
		if (cache) {apply_and_cache_brush(brush, step_sz, num_steps);} else {terrain_hmap_manager_t::apply_brush(brush, step_sz, num_steps);}
		if (cur_tile == NULL) return; // no tile specified, so can't do any updates
		tile_xy_pair const tp(cur_tile->get_tile_xy_pair());

		for (int dy = -1; dy <= 1; ++dy) {
			for (int dx = -1; dx <= 1; ++dx) {
				if (!modified[dy+1][dx+1]) continue;
				tile_t *adj_tile(get_tile_from_xy(tile_xy_pair(tp.x + dx, tp.y + dy)));
				if (adj_tile) {adj_tile->invalidate_mesh_height();}
			}
		}
		cur_tile = NULL;
	}
	virtual bool modify_height_value(int x, int y, hmap_val_t val, bool is_delta, float fract_x, float fract_y) {
		int clamped_x(x), clamped_y(y);
		if (!clamp_xy(clamped_x, clamped_y, fract_x, fract_y)) return 0;
		assert(clamped_x >= 0 && clamped_y >= 0);
		modify_height(tex_mod_map_manager_t::mod_elem_t(clamped_x, clamped_y, val), is_delta); // Note: *not* cached at this level
		if (cur_tile) {cur_tile->fill_adj_mask(modified, x, y);}
		return 1;
	}
};


tiled_terrain_hmap_manager_t terrain_hmap_manager;


bool using_tiled_terrain_hmap_tex() {return (world_mode == WMODE_INF_TERRAIN && terrain_hmap_manager.enabled());}
bool using_hmap_with_detail      () {return (using_tiled_terrain_hmap_tex() && mesh_scale < 0.75);}
float get_tiled_terrain_height_tex(float xval, float yval) {return terrain_hmap_manager.interpolate_height(xval, yval);}
vector3d get_tiled_terrain_height_tex_norm(int x, int y)   {return terrain_hmap_manager.get_norm(x, y);}

bool read_default_hmap_modmap() {

	if (read_hmap_modmap_fn.empty()) return 0;
	if (!terrain_hmap_manager.read_mod(read_hmap_modmap_fn)) return 0;
	cout << "Read heightmap modmap " << read_hmap_modmap_fn << endl;
	return 1;
}

bool write_default_hmap_modmap() {

	if (write_hmap_modmap_fn.empty()) return 0;
	if (!terrain_hmap_manager.write_mod(write_hmap_modmap_fn)) return 0;
	cout << "Wrote heightmap modmap " << write_hmap_modmap_fn << endl;
	return 1;
}


// *** tile_t ***

tile_t::tile_t() : last_occluded_frame(0), weight_tid(0), height_tid(0), shadow_normal_tid(0), size(0), stride(0),
	zvsize(0), gen_tsize(0), decid_trees(tree_data_manager) {}

tile_t::tile_t(unsigned size_, int x, int y) : last_occluded_frame(0), weight_tid(0), height_tid(0), shadow_normal_tid(0),
	size(size_), stride(size+1), zvsize(stride+1), gen_tsize(0), trmax(0.0), min_normal_z(0.0), deltax(DX_VAL), deltay(DY_VAL),
	shadows_invalid(1), recalc_tree_grass_weights(1), mesh_height_invalid(0), in_queue(0), last_occluded(0), has_any_grass(0),
	is_distant(0), no_trees(0), mesh_off(xoff-xoff2, yoff-yoff2), decid_trees(tree_data_manager)
{
	assert(size > 0);
	x1 = x*size;
	y1 = y*size;
	x2 = x1 + size;
	y2 = y1 + size;
	wx1 = x2; wy1 = y2; wx2 = x1; wy2 = y1; // start denormalized
	xstart = get_xval(x1 + mesh_off.dxoff);
	ystart = get_yval(y1 + mesh_off.dyoff);
	radius = calc_radius();
	mzmin  = mzmax = ptzmax = dtzmax = get_camera_pos().z;
	base_tsize = NORM_TEXELS;
}


void tile_t::update_terrain_params() {

	float const off_mult(0.4), height_mult(0.8), dirt_mult(1.0), veg_mult(5.0), off_scale(1.0);
	float const xv1(get_xval(x1)), xv2(xv1 + (x2-x1)*deltax), yv1(get_yval(y1)), yv2(yv1 + (y2-y1)*deltay);

	for (unsigned yp = 0; yp < 2; ++yp) {
		for (unsigned xp = 0; xp < 2; ++xp) {
			terrain_params_t &param(params[yp][xp]);
			float const xv(mesh_scale*(xp ? xv2 : xv1)), yv(mesh_scale*(yp ? yv2 : yv1));
			//param.hoff   = off_scale*eval_mesh_sin_terms(off_mult*xv+123, off_mult*yv+456);
			//param.hscale = min(2.0f, max(0.5f, 0.5f*fabs(eval_mesh_sin_terms(height_mult*xv+789, height_mult*yv+111))));
			float const veg_val(eval_mesh_sin_terms(veg_mult*xv, veg_mult*yv));
			param.veg    = CLIP_TO_01(5.000f*(veg_val + 1.5f));
			param.grass  = CLIP_TO_01(100.0f*(veg_val + 3.0f)); // depends on hoff?
			param.dirt   = CLIP_TO_01(5.0f*(eval_mesh_sin_terms(dirt_mult*xv, dirt_mult*yv) + 1.0f));
		}
	}
}


// used to determine what adjacent tiles modifying this location in global space can affect
void tile_t::fill_adj_mask(bool mask[3][3], int x, int y) const { // mask is {y-1, y, y+1} x {x-1, x, x+1}

	if (x >= x1 && x <= x2 && y >= y1 && y <= y2) {mask[1][1] |= 1;} // ourself
	if (x <= x1) {mask[1][0] |= 1; mask[0][0] |= (y <= y1); mask[2][0] |= (y >= y2);} // left  + top/bottom left  corners
	if (x >= x2) {mask[1][2] |= 1; mask[0][2] |= (y <= y1); mask[2][2] |= (y >= y2);} // right + top/bottom right corners
	mask[0][1] |= (y <= y1); mask[2][1] |= (y >= y2); // top/bottom edges
}


float tile_t::get_min_dist_to_pt(point const &pt, bool xy_only, bool mesh_only) const {

	cube_t const bcube(mesh_only ? get_mesh_bcube() : get_bcube());
	float dsq(0.0);

	for (unsigned i = 0; i < (xy_only ? 2U : 3U); ++i) {
		float const dist(max(0.0f, max((bcube.d[i][0] - pt[i]), (pt[i] - bcube.d[i][1]))));
		dsq += dist*dist;
	}
	return sqrt(dsq);
}

float tile_t::get_max_xy_dist_to_pt(point const &pt) const { // unused

	cube_t const bc(get_mesh_bcube());
	float const dx(max(fabs(pt.x - bc.d[0][0]), fabs(pt.x - bc.d[0][1])));
	float const dy(max(fabs(pt.y - bc.d[1][0]), fabs(pt.y - bc.d[1][1])));
	return sqrt(dx*dx + dy*dy);
}

float tile_t::get_bsphere_radius_inc_water() const {
	return ((is_water_enabled() && mzmax < water_plane_z) ? max(radius, (water_plane_z - 0.5f*(mzmin + mzmax))) : radius); // include water
}


unsigned tile_t::get_gpu_mem() const {

	unsigned mem(pine_trees.get_gpu_mem() + decid_trees.get_gpu_mem() + scenery.get_gpu_mem());
	unsigned const num_texels(stride*stride);
	if (weight_tid        > 0) mem += 4*num_texels; // 4 bytes per texel (RGBA8)
	if (height_tid        > 0) mem += 4*num_texels; // 4 bytes per texel (F32)
	if (shadow_normal_tid > 0) mem += 4*num_texels; // 4 bytes per texel (RGBA8)
	
	for (unsigned i = 0; i < smap_data.size(); ++i) {
		if (smap_data[i].is_allocated()) {mem += 4*shadow_map_sz*shadow_map_sz;} // FBO textures
	}
	return mem;
}


void tile_t::clear() {

	clear_vbo_tid(nullptr);
	tree_map.clear();
	weight_data.clear();
	zvals.clear();
	clear_shadows();
	pine_trees.clear_all();
	decid_trees.clear();
	scenery.clear();
	grass_blocks.clear();
	clear_flowers();
	norm_shadow_data.clear();
}

void tile_t::clear_shadows() {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		smask[l].clear();
		for (unsigned d = 0; d < 2; ++d) {sh_out[l][d].clear();}
	}
	invalidate_shadows();
}

void tile_t::clear_vbo_tid(tile_shadow_map_manager *smap_manager) {

	clear_shadows();
	clear_shadow_map(smap_manager);
	pine_trees.clear_vbos();
	decid_trees.clear_context(); // only necessary if not using instancing
	scenery.clear_vbos();
	flowers.clear_vbo();
	free_texture(weight_tid);
	free_texture(height_tid);
	free_texture(shadow_normal_tid);
	gen_tsize = 0;
}


void setup_height_gen(mesh_xy_grid_cache_t &height_gen, float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values) {

	bool const add_detail(using_hmap_with_detail());

	if (add_detail || !using_tiled_terrain_hmap_tex()) {
		float const xy_scale(add_detail ? HMAP_DETAIL_SCALE : 1.0);
		height_gen.build_arrays(xy_scale*x0, xy_scale*y0, xy_scale*dx, xy_scale*dy, nx, ny, cache_values);
	}
}


void tile_t::create_zvals(mesh_xy_grid_cache_t &height_gen) {

	//RESET_TIME;
	if (enable_terrain_env) {update_terrain_params();}
	zvals.resize(zvsize*zvsize);
	mzmin =  FAR_DISTANCE;
	mzmax = -FAR_DISTANCE;
	float const xy_mult(1.0/float(size)), wpz_max(get_water_z_height() + ocean_wave_height);
	unsigned const block_size(zvsize/4);
	bool const using_hmap(using_tiled_terrain_hmap_tex()), add_detail(using_hmap_with_detail()); // add procedural detail to heightmap
	setup_height_gen(height_gen, get_xval(x1), get_yval(y1), deltax, deltay, zvsize, zvsize, 0); // cache_values=0

	#pragma omp parallel for schedule(static,1)
	for (int y = 0; y < (int)zvsize; ++y) {
		for (unsigned x = 0; x < zvsize; ++x) {
			if (using_hmap) {
				zvals[y*zvsize + x] = terrain_hmap_manager.get_clamped_height((x1 + x), (y1 + y));
				if (add_detail) {zvals[y*zvsize + x] += HMAP_DETAIL_MAG*height_gen.eval_index(x, y, 0);} // less hard-coded - scale by delta between adjacent zvals?
			}
			else {
				float const xv(float(x)*xy_mult), yv(float(y)*xy_mult);
				float const hoff(BILINEAR_INTERP(params, hoff, xv, yv)), hscale(BILINEAR_INTERP(params, hscale, xv, yv));
				zvals[y*zvsize + x] = hoff + hscale*height_gen.eval_index(x, y, 1);
			}
		}
	}
	for (unsigned yy = 0; yy < 4; ++yy) {
		for (unsigned xx = 0; xx < 4; ++xx) {
			sub_zmin[yy][xx] =  FAR_DISTANCE;
			sub_zmax[yy][xx] = -FAR_DISTANCE;
			unsigned const x_end((xx+1)*block_size), y_end((yy+1)*block_size); // last row/column is skipped because it's not rendered
			assert(x_end < zvsize && y_end < zvsize);

			for (unsigned y = yy*block_size; y <= y_end; ++y) {
				for (unsigned x = xx*block_size; x <= x_end; ++x) {
					float const z(zvals[y*zvsize + x]);
					sub_zmin[yy][xx] = min(sub_zmin[yy][xx], z);
					sub_zmax[yy][xx] = max(sub_zmax[yy][xx], z);
				}
			}
		}
	}
	for (vector<float>::const_iterator i = zvals.begin(); i != zvals.end(); ++i) {
		mzmin = min(mzmin, *i);
		mzmax = max(mzmax, *i);
	}
	if (mzmin < wpz_max) { // can have water
		for (unsigned y = 0; y <= size; ++y) {
			for (unsigned x = 0; x <= size; ++x) {
				if (zvals[y*zvsize + x] < wpz_max) { // update water bbox
					wx1 = min(wx1, x1+int(x)); wy1 = min(wy1, y1+int(y));
					wx2 = max(wx2, x1+int(x)); wy2 = max(wy2, y1+int(y));
				}
			}
		}
	}
	assert(mzmin <= mzmax);
	radius = 0.5*sqrt((deltax*deltax + deltay*deltay)*size*size + (mzmax - mzmin)*(mzmax - mzmin));
	ptzmax = dtzmax = mzmin; // no trees yet
	if (!can_have_trees()) {no_trees = 1;} // mark as no_trees so that trees don't pop when water is disabled later
	//PRINT_TIME("Create Zvals");
	if (DEBUG_TILES) {cout << "new tile coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;}
}


// *** shadows + AO lighting ***

void tile_t::calc_mesh_ao_lighting() {

	//RESET_TIME;
	unsigned const NUM_DIRS  = 8; // Note: required to be 8 for adj tile calculation
	unsigned const NUM_STEPS = 8;
	unsigned const ray_len(NUM_STEPS*(NUM_STEPS+1)/2); // 36
	assert(ray_len <= size);

	// caclulate ray step directions and adjacent tiles
	tile_xy_pair const cur_tp(get_tile_xy_pair());
	tile_xy_pair ao_dirs[NUM_DIRS]; // 0  1  2  3  4  5  6  7
	unsigned ix(0);

	for (int y = -1; y <= 1; ++y) {
		for (int x = -1; x <= 1; ++x) {
			if (x != 0 || y != 0) {ao_dirs[ix++] = tile_xy_pair(x, y);}
		}
	}
	assert(ix == NUM_DIRS);

	// create context zvals, which may overlap with other tiles (that need not be created at this point)
	unsigned const context_sz(stride + 2*ray_len);
	vector<float> czv(context_sz*context_sz);
	mesh_xy_grid_cache_t height_gen;
	bool const using_hmap(using_tiled_terrain_hmap_tex()), add_detail(using_hmap_with_detail());
	setup_height_gen(height_gen, get_xval(x1 - ray_len), get_yval(y1 - ray_len), deltax, deltay, context_sz, context_sz, 0); // cache_values=0
	float const dz(0.5*HALF_DXY);
	ao_lighting.resize(stride*stride);

	#pragma omp parallel
	{
		#pragma omp for schedule(static,1)
		for (int y = 0; y < (int)context_sz; ++y) {
			for (int x = 0; x < (int)context_sz; ++x) {
				int const xv(x - ray_len), yv(y - ray_len);
				float &zv(czv[y*context_sz + x]);

				if (xv >= 0 && yv >= 0 && xv < (int)zvsize && yv < (int)zvsize) {
					zv = zvals[yv*zvsize + xv];
				}
				else if (using_hmap) {
					zv = terrain_hmap_manager.get_clamped_height((x1 + xv), (y1 + yv));
					if (add_detail) {zv += HMAP_DETAIL_MAG*height_gen.eval_index(x, y, 0);}
				}
				else {
					zv = height_gen.eval_index(x, y, 1); // Note: not using hoff/hscale here since they are undefined outside the tile bounds
				}
			}
		}

		// calculate ao_lighting values by casting rays through the mesh zvals
		#pragma omp for schedule(static,1)
		for (int y = 0; y < (int)stride; ++y) {
			for (int x = 0; x < (int)stride; ++x) {
				unsigned atten(0);

				for (unsigned d = 0; d < NUM_DIRS; ++d) {
					float z0(zvals[y*zvsize + x]);
					tile_xy_pair step(ao_dirs[d]);
					tile_xy_pair v(x, y);

					for (unsigned s = 0; s < NUM_STEPS; ++s) {
						v    += step;
						z0   += dz;
						//step += step; // multiply by 2 for exponential step size
						step += ao_dirs[d]; // linear increase (Note: must agree with max_ray_length)
						int const xv(v.x + ray_len), yv(v.y + ray_len);
						//assert(xv >= 0 && yv >= 0 && xv < (int)context_sz && yv < (int)context_sz);
						
						if (czv[yv*context_sz + xv] > z0) { // hit a higher point
							atten += (NUM_STEPS - s);
							break;
						}
					} // for s
				} // for d
				assert(atten <= NUM_DIRS*NUM_STEPS);
				float const ao_scale(1.0 - float(atten)/float(NUM_DIRS*NUM_STEPS));
				ao_lighting[y*stride + x] = (unsigned char)(255.0*ao_scale);
			} // for x
		} // for y
	}
	//PRINT_TIME("AO Lighting");
}


void tile_t::calc_shadows_for_light(unsigned l) {

	assert(!smask[l].empty());
	if (is_distant) return; // Note: can be made to work, but won't work as-is
	tile_xy_pair const tp(get_tile_xy_pair());

	// pull from adjacent tiles that already had their shadows calculated
	float const *sh_in[2] = {0, 0};
	point const lpos(get_light_pos(l));
	tile_xy_pair const adj_tp[2] = {tile_xy_pair((tp.x + ((lpos.x < 0.0) ? -1 : 1)), tp.y),
									tile_xy_pair(tp.x, (tp.y + ((lpos.y < 0.0) ? -1 : 1)))}; // toward the light source

	for (unsigned d = 0; d < 2; ++d) { // d = tile adjacency dimension, shared edge is in !d
		sh_out[l][!d].resize(zvsize, MESH_MIN_Z); // init value really should not be used, but it sometimes is
		tile_t *adj_tile(get_tile_from_xy(adj_tp[d]));
		if (adj_tile == NULL || adj_tile->is_distant) continue; // no adjacent tile
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


void tile_t::proc_tile_queue(tile_t *init_tile, unsigned l) {

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
			if (adj_tile == NULL || adj_tile->is_distant || adj_tile->smask[l].empty() || adj_tile->in_queue) continue; // no adjacent tile, not initialized, or already in queue
			tile_queue.push_front(adj_tile); // changed, push to adjacent tiles
			adj_tile->in_queue = 1;
		}
	}
}


void tile_t::calc_shadows(bool calc_sun, bool calc_moon, bool no_push) {

	bool calc_light[NUM_LIGHT_SRC] = {0};
	calc_light[LIGHT_SUN ] = calc_sun;
	calc_light[LIGHT_MOON] = calc_moon;

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) { // calculate mesh shadows for each light source
		if (!calc_light[l])    continue; // light not enabled
		if (!smask[l].empty()) continue; // already calculated (cached)
		smask[l].resize(zvals.size(), 0);
		//if (normal_zmin < 1.0 && get_light_pos(l).get_norm().xy_mag() < normal_zmin) { // terrain slope lower than sun slope
		if (no_push) {calc_shadows_for_light(l);} else {proc_tile_queue(this, l);}
	}
}


void tile_t::push_tree_ao_shadow(int dx, int dy, point const &pos, float tradius) const {

	tile_t *const adj_tile(get_adj_tile_smap(dx, dy));
	if (!adj_tile || adj_tile->is_distant) return;
	point const pos2(pos + mesh_off.subtract_from(adj_tile->mesh_off));
	adj_tile->add_tree_ao_shadow(pos2, tradius, 1);
}


void tile_t::add_tree_ao_shadow(point const &pos, float tradius, bool no_adj_test) {

	int const xc(round_fp((pos.x - xstart)/deltax)), yc(round_fp((pos.y - ystart)/deltay));
	int rval(max(int(tradius/deltax), int(tradius/deltay)) + 1);
	int const x1(max(0, xc-rval)), y1(max(0, yc-rval)), x2(min((int)size, xc+rval)), y2(min((int)size, yc+rval));
	float const scale(0.6/rval);
	bool updated(0);

	for (int y = y1; y <= y2; ++y) {
		for (int x = x1; x <= x2; ++x) {
			float const dx(abs(x - xc)), dy(abs(y - yc)), dist(sqrt(dx*dx + dy*dy));
			if (dist < rval) {tree_map[y*stride + x] *= (0.2 + 0.8*scale*dist); updated = 1;}
		}
	}
	if (!no_adj_test) {
		bool const x_test[3] = {(xc <= rval), 1, (xc >= (int)size-rval)};
		bool const y_test[3] = {(yc <= rval), 1, (yc >= (int)size-rval)};

		for (int dy = -1; dy <= 1; ++dy) {
			for (int dx = -1; dx <= 1; ++dx) {
				if (dx == 0 && dy == 0) continue;
				if (x_test[dx+1] && y_test[dy+1]) {push_tree_ao_shadow(dx, dy, pos, tradius);}
			}
		}
	}
	if (updated) { // tree_map was modified, so we need to recalculate both weights (grass replaced with dirt) and shadows
		invalidate_shadows();
		if (has_any_grass) {recalc_tree_grass_weights = 1;} // Note: may be slow, and doesn't have a big impact
	}
}


void tile_t::apply_ao_shadows_for_trees(tile_t const *const tile, bool no_adj_test) {

	if (pine_trees_enabled() && can_have_trees()) {
		assert(pine_trees_generated());
		point const pt_off(tile->ptree_off.subtract_from(mesh_off));

		for (small_tree_group::const_iterator i = tile->pine_trees.begin(); i != tile->pine_trees.end(); ++i) {
			add_tree_ao_shadow((i->get_pos() + pt_off), i->get_pine_tree_radius(), no_adj_test);
		}
	}
	if (decid_trees_enabled() && can_have_trees()) {
		assert(decid_trees.was_generated());
		point const dt_off(tile->dtree_off.subtract_from(mesh_off));

		for (tree_cont_t::const_iterator i = tile->decid_trees.begin(); i != tile->decid_trees.end(); ++i) {
			add_tree_ao_shadow((i->get_center() + dt_off), 0.6*i->get_radius(), no_adj_test); // less dense => smaller radius
		}
	}
	if (!no_adj_test && !is_distant) { // pull mode
		for (int dy = -1; dy <= 1; ++dy) {
			for (int dx = -1; dx <= 1; ++dx) {
				if (dx == 0 && dy == 0) continue;
				tile_t const *const adj_tile(get_adj_tile_smap(dx, dy));
				if (adj_tile && !adj_tile->is_distant) {apply_ao_shadows_for_trees(adj_tile, 1);}
			}
		}
	}
}


void tile_t::apply_tree_ao_shadows() { // should this generate a float or unsigned char shadow weight instead?

	if (is_distant) return; // not needed/used
	tree_map.resize(0);
	tree_map.resize(stride*stride, 255);
	bool const no_adj_test(trmax < min(deltax, deltay));
	apply_ao_shadows_for_trees(this, no_adj_test);
}


void tile_t::check_shadow_map_and_normal_texture() {

	if (shadow_normal_tid && !shadows_invalid) return; // up-to-date
	//RESET_TIME;
	bool const tid_is_valid(shadow_normal_tid != 0);
	if (!tid_is_valid) {setup_texture(shadow_normal_tid, 0, 0, 0, 0, 0);}
	bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6), mesh_shadows(mesh_shadows_enabled());
	assert(has_sun || has_moon);
	if (mesh_shadows) {calc_shadows(has_sun, has_moon);}
	//PRINT_TIME("Calc Shadows");
	if (enable_tiled_mesh_ao && ao_lighting.empty()) {calc_mesh_ao_lighting();}
	upload_shadow_map_and_normal_texture(tid_is_valid);
	shadows_invalid = 0;
	//PRINT_TIME("Calc and Upload Shadows + AO");
}


void tile_t::upload_shadow_map_and_normal_texture(bool tid_is_valid) {

	bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6), mesh_shadows(mesh_shadows_enabled());
	vector<norm_comp_with_shadow> &data(norm_shadow_data); // stored as (n.x, n.y, ao_lighting, shadow_val)
	bool const init_data(data.empty());

	if (init_data) {
		data.resize(stride*stride);
		min_normal_z = 1.0;
	}
	for (unsigned y = 0; y < stride; ++y) {
		for (unsigned x = 0; x < stride; ++x) {
			unsigned const ix(y*stride + x), ix2(y*zvsize + x);

			if (init_data) {
				vector3d const norm(get_norm(ix2));
				min_normal_z = min(min_normal_z, norm.z);
				UNROLL_2X(data[ix].v[i_] = (unsigned char)(127.0*(norm[i_] + 1.0));); // Note: we only set x and y here, z is calculated in the shader
			}
			unsigned char shadow_val(tree_map.empty() ? 255 : tree_map[ix]); // fully lit (if not nearby trees)
			// 67% ambient if AO lighting is disabled (to cancel out with the scale by 1.5 in the shaders)
			data[ix].v[2] = (ao_lighting.empty() ? 170 : ao_lighting[ix]); // Note: must always do this so that adj tiles can update tree AO at tile borders
			data[ix].v[2] = (unsigned char)(data[ix].v[2] * (0.3 + 0.7*shadow_val/255.0)); // add ambient occlusion from trees
			shadow_val = 128 + shadow_val/2; // divide tree shadow effect by 2x to offset for shadow map effect

			if (!mesh_shadows) {
				// do nothing
			}
			else if (has_sun && has_moon) {
				bool const no_sun( (smask[LIGHT_SUN ][ix2] & SHADOWED_ALL) != 0);
				bool const no_moon((smask[LIGHT_MOON][ix2] & SHADOWED_ALL) != 0);
				shadow_val *= blend_light(light_factor, !no_sun, !no_moon);
			}
			else if (smask[has_sun ? LIGHT_SUN : LIGHT_MOON][ix2] & SHADOWED_ALL) {
				shadow_val = 0; // fully in shadow
			}
			data[ix].v[3] = shadow_val;
		}
	}
	bind_2d_texture(shadow_normal_tid);

	if (tid_is_valid) { // overwrite old data
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, stride, stride, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
	}
	else { // allocate and write
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, stride, stride, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
	}
}


tile_smap_data_t tile_shadow_map_manager::new_smap_data(unsigned tu_id, tile_t *tile, unsigned light) {
	assert(tile != nullptr);
	assert(light < NUM_LIGHT_SRC);
	if (free_list[light].empty()) {return tile_smap_data_t(tu_id, shadow_map_sz, tile);}
	smap_data_state_t const state(free_list[light].back());
	free_list[light].pop_back();
	return tile_smap_data_t(tu_id, shadow_map_sz, tile, state);
}

void tile_shadow_map_manager::release_smap_data(tile_smap_data_t &smd, unsigned light) {
	if (!smd.is_allocated()) return;
	assert(light < NUM_LIGHT_SRC);
	free_list[light].push_back(smd); // adds the base class, which contains the GL state
	smd.disown(); // no longer owned by smd
}

void tile_shadow_map_manager::clear_context() {
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		for (auto i = free_list[l].begin(); i != free_list[l].end(); ++i) {i->free_gl_state();}
		free_list[l].clear();
	}
}


void tile_t::setup_shadow_maps(tile_shadow_map_manager &smap_manager) {

	if (!shadow_map_enabled()) return; // disabled

	if (get_dist_to_camera_in_tiles(1) < SMAP_NEW_THRESH && smap_data.empty()) { // allocate new shadow maps
		for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) { // uses tu_id 13 and 14
			smap_data.push_back(smap_manager.new_smap_data(13+i, this, i));
		}
	}
	cube_t bcube(get_bcube());
	// extend bcube upwards to include any models above the mesh that cast shadows on this tile
	// FIXME: still not correct for low sun pos - need a more accurate way to determine which models can shadow this tile
	cube_t const models_bcube(get_all_models_bcube());
	if (models_bcube != all_zeros_cube) {bcube.d[2][1] = max(bcube.d[2][1], models_bcube.d[2][1]);}
	smap_data.create_if_needed(bcube);
}

void tile_t::clear_shadow_map(tile_shadow_map_manager *smap_manager) {
	
	if (smap_data.empty()) return;
	assert(smap_data.size() == NUM_LIGHT_SRC);

	if (smap_manager) {
		for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) {smap_manager->release_smap_data(smap_data[i], i);}
	}
	smap_data.clear(); // also frees any leftover GL state
}


// *** mesh creation ***


void tile_t::ensure_height_tid() {

	if (height_tid || is_distant) return; // already exists, or tile is distant and height_tid is unnecessary
	assert(zvals.size() == zvsize*zvsize);
	vector<float> data(stride*stride);

	for (unsigned y = 0; y < stride; ++y) {
		for (unsigned x = 0; x < stride; ++x) {
			data[y*stride+x] = zvals[y*zvsize+x]; // remove the last column
		}
	}
	setup_texture(height_tid, 0, 0, 0, 0, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, stride, stride, 0, GL_RED, GL_FLOAT, &data.front());
}


void tile_t::create_texture(mesh_xy_grid_cache_t &height_gen) {

	assert(zvals.size() == zvsize*zvsize);
	//RESET_TIME;
	unsigned const tsize(stride);
	int sand_tex_ix(-1), dirt_tex_ix(-1), grass_tex_ix(-1), rock_tex_ix(-1), snow_tex_ix(-1);

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

	if (weight_tid == 0) { // create weights
		has_any_grass = 0;
		grass_blocks.clear();
		weight_data.resize(4*tsize*tsize); // RGBA
		unsigned const grass_block_dim(get_grass_block_dim());
		float const xy_mult(1.0/float(size)), water_level(get_water_z_height());
		float const MESH_NOISE_SCALE = 0.003;
		float const MESH_NOISE_FREQ  = 80.0;
		float const dz_inv(1.0/(zmax - zmin));
		float const noise_scale(((mesh_gen_shape == 2) ? 2.0 : 1.0)*MESH_NOISE_SCALE*mesh_scale_z); // add more noise for ridged
		float const steep_mult_grass(1.0/(sthresh[0][1] - sthresh[0][0]));
		float const steep_mult_snow (1.0/(sthresh[1][1] - sthresh[1][0]));
		float const steep_mult_rock (1.0/(0.8f*sthresh[0][0] - 0.5f*sthresh[0][0]));
		int k1, k2, k3, k4;
		height_gen.build_arrays(MESH_NOISE_FREQ*get_xval(x1), MESH_NOISE_FREQ*get_yval(y1), MESH_NOISE_FREQ*deltax,
			MESH_NOISE_FREQ*deltay, tsize, tsize, 0, 1); // force_sine_mode=1

		for (unsigned y = 0; y < tsize-DEBUG_TILE_BOUNDS; ++y) { // not threadsafe
			float const yv(float(y)*xy_mult);

			for (unsigned x = 0; x < tsize-DEBUG_TILE_BOUNDS; ++x) {
				float weights[NTEX_DIRT] = {0};
				unsigned const ix(y*zvsize + x);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const rand_offset(noise_scale*height_gen.eval_index(x, y, 0, 50));
				float const mhmin(min(min(mh00, mh01), min(mh10, mh11))), mhmax(max(max(mh00, mh01), max(mh10, mh11)));
				float const relh1(relh_adj_tex + (mhmin - zmin)*dz_inv + rand_offset);
				float const relh2(relh_adj_tex + (mhmax - zmin)*dz_inv + rand_offset);
				get_tids(relh1, k1, k2);
				get_tids(relh2, k3, k4);
				bool const same_tid(k1 == k4);
				float t(0.0);
				k2 = k4;
			
				if (!same_tid) {
					float const relh(relh_adj_tex + (mh00 - zmin)*dz_inv);
					get_tids(relh, k1, k2, &t);
				}
				float weight_scale(1.0);
				bool const grass(lttex_dirt[k1].id == GROUND_TEX || lttex_dirt[k2].id == GROUND_TEX), snow(lttex_dirt[k2].id == SNOW_TEX);
				has_any_grass |= grass;

				if (grass || snow) {
					float const *const sti(sthresh[snow]);
					vector3d const normal(get_norm_not_normalized(ix));
					float const vnz(normal.z/normal.mag());

					if (vnz < sti[1]) { // handle steep slopes (dirt/rock texture replaces grass texture)
						if (grass) { // ground/grass
							float rock_weight((lttex_dirt[k1].id == GROUND_TEX || lttex_dirt[k2].id == ROCK_TEX) ? t : 0.0);
							float const steepness(1.0 - CLIP_TO_01((vnz - 0.5f*sti[0])*steep_mult_rock));
							rock_weight  = rock_weight*(1.0 - steepness) + steepness;
							weight_scale = CLIP_TO_01((vnz - sti[0])*steep_mult_grass);
							weights[rock_tex_ix] += (1.0 - weight_scale)*rock_weight;
							weights[dirt_tex_ix] += (1.0 - weight_scale)*(1.0 - rock_weight);
						}
						else { // snow
							weight_scale = CLIP_TO_01(2.0f*(vnz - sti[0])*steep_mult_snow);
							weights[rock_tex_ix] += 1.0 - weight_scale;
						}
					}
				}
				weights[k2] += weight_scale*t;
				weights[k1] += weight_scale*(1.0 - t);
				unsigned const ix_val(y*tsize + x), off(4*ix_val);
				float const xv(float(x)*xy_mult);
				float const dirt_scale(BILINEAR_INTERP(params, dirt, xv, yv));

				if (dirt_scale < 1.0) { // apply dirt scale: convert dirt to sand
					weights[sand_tex_ix ] += (1.0 - dirt_scale)*weights[dirt_tex_ix];
					weights[dirt_tex_ix ] *= dirt_scale;
				}
				if (grass) {
					float const grass_scale((mhmin < water_level) ? 0.0 : BILINEAR_INTERP(params, grass, xv, yv)); // no grass under water

					if (grass_scale < 1.0) { // apply grass scale: convert grass to sand
						float const gscale(CLIP_TO_01(2.5f*(grass_scale - 0.5f) + 0.5f));
						weights[sand_tex_ix ] += (1.0 - gscale)*weights[grass_tex_ix];
						weights[grass_tex_ix] *= gscale;
					}
					if (!is_distant && x < size && y < size && gen_grass_map()) {
						unsigned const bx(x/GRASS_BLOCK_SZ), by(y/GRASS_BLOCK_SZ), bix(by*grass_block_dim + bx);
						if (grass_blocks.empty()) {grass_blocks.resize(grass_block_dim*grass_block_dim);}
						assert(bix < grass_blocks.size());
						grass_block_t &gb(grass_blocks[bix]);
					
						if (gb.ix == 0) { // not yet set
							gb.ix = (((x1 + x) + 1567*(y1 + y)) % NUM_RND_GRASS_BLOCKS) + 1; // select a random block
							//gb.ix = (rand() % NUM_RND_GRASS_BLOCKS) + 1; // select a random block
							gb.zmin = mhmin;
							gb.zmax = mhmax;
						}
						else {
							gb.zmin = min(gb.zmin, mhmin);
							gb.zmax = max(gb.zmax, mhmax);
						}
					}
				} // end grass
				for (unsigned i = 0; i < NTEX_DIRT-1; ++i) { // Note: weights should sum to 1.0, so we can calculate w4 as 1.0-w0-w1-w2-w3
					weight_data[off+i] = (weights[i] ? (unsigned char)(255.0*CLIP_TO_01(weights[i])) : 0);
				}
			} // for x
		} // for y
	}
	else { // use existing weights
		assert(recalc_tree_grass_weights); // can only get here in this case
		assert(weight_data.size() == 4*tsize*tsize);
	}
	vector<unsigned char> data(weight_data); // deep copy so that tree_map doesn't alter original weights

	if (!tree_map.empty()) {
		for (unsigned i = 0; i < tsize*tsize; ++i) { // replace grass under trees with dirt
			unsigned const off(4*i);
			if (tree_map[i] == 255 || data[off+grass_tex_ix] == 0) continue; // no trees or no grass
			float const v(tree_map[i]/255.0);
			data[off+dirt_tex_ix]   = (unsigned char)(max(0.0, min(255.0, (data[off+dirt_tex_ix] + (1.0 - v)*data[off+grass_tex_ix]))));
			data[off+grass_tex_ix] *= v;
		}
	}
	if (weight_tid == 0) { // create weight texture
		setup_texture(weight_tid, 0, 0, 0, 0, 0);
		assert(weight_tid > 0 && glIsTexture(weight_tid));
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tsize, tsize, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data.front()); // internal_format = GL_COMPRESSED_RGBA - too slow
	}
	else { // update texture
		bind_2d_texture(weight_tid);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tsize, tsize, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
	}
	recalc_tree_grass_weights = 0;
	//PRINT_TIME("Texture Upload");
}


bool tile_t::update_range(tile_shadow_map_manager &smap_manager) { // if returns 0, tile will be deleted

	update_pine_tree_state(0); // can free pine tree vbos
	float const dist(get_rel_dist_to_camera());
	if (dist > CLEAR_DIST_TILES || mesh_height_invalid) {clear_vbo_tid(&smap_manager);}
	if (dist*TILE_RADIUS > SMAP_DEL_THRESH) {clear_shadow_map(&smap_manager);} // too far, delete old shadow maps
	return (dist < DELETE_DIST_TILES && !mesh_height_invalid);
}


// *** trees ***

void tile_t::init_pine_tree_draw() {

	if (!can_have_trees()) return; // no pine trees (yet)
	float const density[4] = {params[0][0].veg, params[0][1].veg, params[1][0].veg, params[1][1].veg};
	ptree_off.set_from_xyoff2();
	if (enable_instanced_pine_trees()) {pine_trees.instanced = 1;}
	pine_trees.gen_trees(x1+ptree_off.dxoff, y1+ptree_off.dyoff, x2+ptree_off.dxoff, y2+ptree_off.dyoff, density);
	pine_trees.calc_trunk_pts();
	postproc_trees(pine_trees, ptzmax);
}


void tile_t::update_pine_tree_state(bool upload_if_needed) {

	if (!pine_trees_enabled() || pine_trees.empty()) return;
	float const weight(get_tree_far_weight());
	float const weights[2] = {1.0-weight, weight}; // {high, low} detail

	for (unsigned d = 0; d < 2; ++d) {
		if (weights[d] > 0.0) { // needed
			if (upload_if_needed) {pine_trees.finalize_upload_and_clear_pts(d != 0);} // needed for drawing
		}
		else { // not needed
			pine_trees.clear_vbo_and_ids_if_needed(d != 0);
		}
	}
}


// non-const because of cached instance data within pine_trees
void tile_t::draw_tree_leaves_lod(vector3d const &xlate, bool low_detail, int xlate_loc) {

	bool const draw_all(low_detail || camera_pdu.point_visible_test(get_center())); // tile center is in view
	pine_trees.draw_pine_leaves(0, low_detail, draw_all, contains_camera(), xlate, xlate_loc);
}


// Note: xlate has different meanings here: for near leaves, it's a vec3 attribute; for far leaves, it's a vec2 uniform; for branches, it's -1
void tile_t::draw_pine_trees(shader_t &s, vector<vert_wrap_t> &trunk_pts, bool draw_branches, bool draw_near_leaves,
	bool draw_far_leaves, bool reflection_pass, int xlate_loc)
{
	if (pine_trees.empty() || !can_have_trees()) return;
	point const camera(get_camera_pos());
	vector3d const xlate(ptree_off.get_xlate());
	vector2d const camera_xlate(xlate.x-camera.x, xlate.y-camera.y);
	fgPushMatrix();
	translate_to(xlate);
	
	if (draw_branches) {
		float const dscale(get_tree_dist_scale());

		if (dscale < 1.0) { // close, draw as polygons
			pine_trees.draw_branches(0, xlate, &trunk_pts);
		}
		else if (dscale < 2.0 && get_tree_far_weight() < 0.5) { // far away, use low detail branches
			pine_trees.add_trunk_pts(xlate, trunk_pts);
		} // else very far, skip branches
	}
	if (draw_near_leaves || draw_far_leaves) { // could use reflection_pass as an optimization
		float const weight(1.0 - get_tree_far_weight()); // 0 => low detail, 1 => high detail
		
		if (draw_far_leaves && weight < 1.0) {
			assert(xlate_loc >= 0);
			s.set_uniform_vector2d(xlate_loc, camera_xlate);
		}
		if (weight > 0.0 && weight < 1.0) { // use geomorphing with dithering (since alpha doesn't blend in the correct order)
			if (draw_near_leaves) {
				int const loc(s.get_uniform_loc("max_noise"));
				s.set_uniform_float(loc, weight);
				draw_tree_leaves_lod(xlate, 0, xlate_loc); // near leaves
				s.set_uniform_float(loc, 1.0);
			}
			if (draw_far_leaves) {
				int const loc(s.get_uniform_loc("min_noise"));
				s.set_uniform_float(loc, weight);
				draw_tree_leaves_lod(xlate, 1, xlate_loc); // far leaves
				s.set_uniform_float(loc, 0.0);
			}
		}
		else if ((weight == 0.0) ? draw_far_leaves : draw_near_leaves) {
			draw_tree_leaves_lod(xlate, (weight == 0.0), xlate_loc);
		}
	}
	fgPopMatrix();
}


void tile_t::gen_decid_trees_if_needed() {

	if (num_trees > 0 && max_unique_trees == 0) {
		cout << "Warning: max_unique_trees needs to be set to something reasonable for tiled terrain mode trees to work efficiently. Setting to 100." << endl;
		max_unique_trees = 100;
	}
	if (decid_trees.was_generated() || !can_have_trees()) return; // already generated, or distant tile (no trees yet)
	assert(decid_trees.empty());
	dtree_off.set_from_xyoff2();
	decid_trees.gen_deterministic(x1+dtree_off.dxoff, y1+dtree_off.dyoff, x2+dtree_off.dxoff, y2+dtree_off.dyoff, vegetation*get_avg_veg());
	postproc_trees(decid_trees, dtzmax);
}


void tile_t::update_decid_trees() {
	if (decid_trees_enabled()) {decid_trees.check_render_textures();}
}


void tile_t::draw_decid_trees(shader_t &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves, bool reflection_pass) {

	if (decid_trees.empty() || !can_have_trees()) return;
	// Note: shadow_only is always zero here since it's not really correct and doesn't help performance
	decid_trees.draw_branches_and_leaves(s, lod_renderer, draw_branches, draw_leaves, 0, reflection_pass, dtree_off.get_xlate());
}


// *** scenery/grass ***

void tile_t::update_scenery() {

	if (!scenery_enabled() || is_distant) return; // no scenery
	float const dist_scale(get_scenery_dist_scale(0)); // tree_dist_scale should correlate with mesh scale
	if (scenery.generated && dist_scale > 1.2) {scenery.clear();} // too far away
	if (scenery.generated || dist_scale > 1.0 || !is_visible()) return; // already generated, too far away, or not visible
	scenery_off.set_from_xyoff2();
	scenery.gen(x1+scenery_off.dxoff, y1+scenery_off.dyoff, x2+scenery_off.dxoff, y2+scenery_off.dyoff, vegetation*get_avg_veg(), 1);
}


void tile_t::draw_scenery(shader_t &s, bool draw_opaque, bool draw_leaves, bool reflection_pass) {

	if (!scenery.generated || get_scenery_dist_scale(reflection_pass) > 1.0) return;
	fgPushMatrix();
	vector3d const xlate(scenery_off.get_xlate());
	translate_to(xlate);
	float const scale_val(get_scenery_thresh(reflection_pass)*get_tile_width());
	if (draw_opaque) {scenery.draw_opaque_objects(s, 0, xlate, 0, scale_val, reflection_pass);} // shader not passed in here
	if (draw_leaves) {scenery.draw_plant_leaves  (s, 0, xlate, reflection_pass);}
	fgPopMatrix();
}


void tile_t::pre_draw_grass_flowers(shader_t &s, bool use_cloud_shadows) const {

	bind_texture_tu(height_tid, 2);
	bind_texture_tu(shadow_normal_tid, 4);
	shader_shadow_map_setup(s);
	if (use_cloud_shadows) {s.add_uniform_vector3d("cloud_offset", vector3d(get_xval(x1), get_yval(y1), 0.0));}
	s.add_uniform_vector2d("xlate", vector2d(get_xval(x1 + xoff - xoff2), get_yval(y1 + yoff - yoff2)));
}


void tile_t::draw_grass(shader_t &s, vector<vector<vector2d> > *insts, bool use_cloud_shadows, int lt_loc) {

	if (grass_blocks.empty()) return; // or can test has_any_grass
	float const grass_thresh(get_grass_thresh_pad());
	point const camera(get_camera_pos());
	if (get_min_dist_to_pt(camera) > grass_thresh) return; // too far away to draw
	pre_draw_grass_flowers(s, use_cloud_shadows);
	bind_texture_tu(weight_tid, 3);
	unsigned const grass_block_dim(get_grass_block_dim());
	assert(grass_blocks.size() == grass_block_dim*grass_block_dim);
	float const llcx(get_xval(x1 + xoff - xoff2)), llcy(get_yval(y1 + yoff - yoff2));
	float const dx_step(GRASS_BLOCK_SZ*deltax), dy_step(GRASS_BLOCK_SZ*deltay);
	float const lod_scale(GRASS_LOD_SCALE/get_scaled_tile_radius());
	float const block_grass_thresh(grass_thresh + (SQRT2*radius)/grass_block_dim);
	point const adj_camera(camera + point(0.0, 0.0, 2.0*grass_length));

	for (unsigned y = 0; y < grass_block_dim; ++y) {
		for (unsigned x = 0; x < grass_block_dim; ++x) {
			grass_block_t const &gb(grass_blocks[y*grass_block_dim+x]);
			if (gb.ix == 0) continue; // empty block
			cube_t const bcube(llcx+x*dx_step, llcx+(x+1)*dx_step, llcy+y*dy_step, llcy+(y+1)*dy_step, gb.zmin, (gb.zmax + grass_length));
			point const close_pt(bcube.closest_pt(camera));
			if (!dist_less_than(camera, close_pt, block_grass_thresh) || !camera_pdu.cube_visible(bcube)) continue;
			bool back_facing(1);

			for (unsigned yy = y*GRASS_BLOCK_SZ; yy <= (y+1)*GRASS_BLOCK_SZ && back_facing; ++yy) {
				for (unsigned xx = x*GRASS_BLOCK_SZ; xx <= (x+1)*GRASS_BLOCK_SZ && back_facing; ++xx) {
					unsigned const ix(yy*zvsize + xx);
					back_facing &= (dot_product(get_norm_not_normalized(ix), (adj_camera - point(llcx+xx*deltax, llcy+yy*deltay, zvals[ix]))) < 0.0);
				}
			}
			if (back_facing) continue;
			unsigned const lod_level(min(NUM_GRASS_LODS-1, unsigned(lod_scale*distance_to_camera(close_pt))));
			assert(insts != NULL);
			insts[lod_level].resize(NUM_RND_GRASS_BLOCKS); // may already be the correct size
			unsigned const bix(gb.ix - 1);
			assert(bix < NUM_RND_GRASS_BLOCKS);
			insts[lod_level][bix].push_back(vector2d(x*dx_step, y*dy_step));
		} // for x
	} // for y
	for (unsigned lod = 0; lod < NUM_GRASS_LODS; ++lod) {
		for (unsigned bix = 0; bix < insts[lod].size(); ++bix) {
			vector<vector2d> &v(insts[lod][bix]);
			if (v.empty()) continue;
			glVertexAttribPointer(lt_loc, 2, GL_FLOAT, GL_FALSE, sizeof(vector2d), get_dynamic_vbo_ptr(&v.front(), v.size()*sizeof(vector2d)));
			grass_tile_manager.render_block(bix, lod, 1.0, v.size());
			v.clear();
		} // for bix
	} // for lod
}


void tile_t::draw_flowers(shader_t &s, bool use_cloud_shadows) {

	if (grass_blocks.empty()) return; // no grass, no flowers
	float const flower_thresh(FLOWER_REL_DIST*get_grass_thresh_pad());
	if (get_min_dist_to_pt(get_camera_pos()) > flower_thresh) return; // too far away to draw
	flowers.gen_flowers(weight_data, stride, x1-xoff2, y1-yoff2);
	if (flowers.empty()) return; // no flowers generated
	pre_draw_grass_flowers(s, use_cloud_shadows);
	flowers.check_vbo();
	flowers.draw_triangles(s);
}


// *** rendering ***


void tile_t::pre_draw(mesh_xy_grid_cache_t &height_gen) {

	if (tree_map.empty() && any_trees_enabled()) {apply_tree_ao_shadows();}
	if (weight_tid == 0 || recalc_tree_grass_weights) {create_texture(height_gen);}
	check_shadow_map_and_normal_texture();
	ensure_height_tid();
}


unsigned tile_t::get_lod_level(bool reflection_pass) const {

	unsigned lod_level((reflection_pass && min_normal_z > 0.1) ? min(NUM_LODS-1, 1U) : 0);
	float dist(get_dist_to_camera_in_tiles(0));
	
	if (min_normal_z > 0.0) { // normals have been calculated, adjust detail based on max slope
		if (min_normal_z > 0.9) { // flat, lower detail
			if (!is_water_enabled() || mzmax < water_plane_z || mzmin > water_plane_z) { // use get_water_z_height()?
				dist *= (min_normal_z > 0.95) ? 4.0 : 2.0;
			}
		}
		else if (min_normal_z < 0.25) { // high slope, higher detail
			dist /= -log(2.0*min_normal_z)/log(2.0);
		}
	}
    while (dist > (reflection_pass ? 1.8 : 2.0) && lod_level+1 < NUM_LODS && (size>>(lod_level+1)) >= 4) { // never divide smaller than 4x4 quads
        dist /= 2;
        ++lod_level;
    }
	assert(lod_level < NUM_LODS);
	return lod_level;
}


void tile_t::shader_shadow_map_setup(shader_t &s, xform_matrix const *const mvm) const {

	// Note: some part of this call is shared across all tiles; however, in the case where more than one smap light is enabled,
	// the tu_id and enables may alternate between values for each tile, requiring every uniform to be reset per tile anyway
	smap_data.set_for_all_lights(s, mvm);
}


void tile_t::bind_textures() const {

	assert(weight_tid > 0 && height_tid > 0 && shadow_normal_tid > 0);
	bind_2d_texture(weight_tid);
	bind_texture_tu(height_tid, 12);
	bind_texture_tu(shadow_normal_tid, 7);
}


void crack_ibuf_t::gen_offsets(vector<unsigned> &indices, unsigned size) {

	offsets.clear();
	offsets.resize(4*NUM_LODS*NUM_LODS);

	for (unsigned dim = 0; dim < 2; ++dim) { // x,y
		for (unsigned dir = 0; dir < 2; ++dir) { // lo, hi
			for (unsigned cur_lod = 0; cur_lod < NUM_LODS; ++cur_lod) {
				for (unsigned adj_lod = 0; adj_lod < NUM_LODS; ++adj_lod) {
					ix_sz_pair &ixsz(offsets[get_index(dim, dir, cur_lod, adj_lod)]);
					ixsz.ix = indices.size();
					unsigned const step(1 << cur_lod), stride(size+1);

					for (unsigned lod = adj_lod; lod > cur_lod; --lod) { // for all levels of high=>low LOD transitions
						unsigned const lo_step(1 << lod), hi_step(1 << (lod - 1));
						unsigned const size_lo_ceil((size + step - 1)/lo_step);

						for (unsigned xy = 0, nxy = 0; nxy < size_lo_ceil; xy += lo_step, ++nxy) {
							for (unsigned n = 0; n < 3; ++n) { // one triangle
								if (dim == 0) { // adjacent in x, step in y
									indices.push_back(min((xy + n*hi_step), size)*stride + dir*size);
								}
								else { // adjacent in y, step in x
									indices.push_back(dir*size*stride + min((xy + n*hi_step), size));
								}
							} // for n
						} // for xy
					} // for lod
					ixsz.sz = indices.size() - ixsz.ix; // often 0
				} // for adj_lod
			} // for cur_lod
		} // for dir
	} // for dim
}


unsigned crack_ibuf_t::get_index(unsigned dim, unsigned dir, unsigned cur_lod, unsigned adj_lod) const {
	assert(dim < 2 && dir < 2 && cur_lod < NUM_LODS && adj_lod < NUM_LODS);
	return (NUM_LODS*(NUM_LODS*(2*dim + dir) + cur_lod) + adj_lod);
}


void tile_t::draw(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], crack_ibuf_t const &crack_ibuf, bool reflection_pass) const {

	assert(size > 0);
	fgPushMatrix();
	vector3d const xlate(mesh_off.get_xlate() + vector3d(xstart, ystart, 0.0));
	translate_to(xlate); // Note: not easy to replace with a uniform, due to texgen and fog dist calculations in the shader
	bind_textures();
	bool const draw_near_water(!is_distant && !reflection_pass && is_water_enabled() && has_water());

	if (!reflection_pass && cloud_shadows_enabled()) {
		s.add_uniform_vector3d("cloud_offset", vector3d(get_xval(x1), get_yval(y1), 0.0));
	}
	if (draw_near_water) {s.add_uniform_vector3d("tc_xlate", xlate);} // for underwater caustics texture
	shader_shadow_map_setup(s);
	unsigned const lod_level(get_lod_level(reflection_pass)), num_ixs(ivbo_ixs[lod_level+1] - ivbo_ixs[lod_level]);
	vbo_mgr.pre_render();
	vert_wrap_t::set_vbo_arrays(0); // normals are stored in shadow_normal_tid, tex coords come from texgen, color is constant
	glDrawRangeElements(GL_TRIANGLE_STRIP, 0, stride*stride, num_ixs, GL_UNSIGNED_INT, (void *)(ivbo_ixs[lod_level]*sizeof(unsigned)));
	
	// draw cracks
	for (unsigned dim = 0; dim < 2; ++dim) { // x,y
		for (unsigned dir = 0; dir < 2; ++dir) { // lo, hi
			int dx(0), dy(0);
			(dim ? dy : dx) += (dir ? 1 : -1);
			tile_t *adj(get_adj_tile(dx, dy)); // FIXME: handle (is_distant != adj->is_distant)
			if (adj == NULL) continue; // no adjacent tile
			//if (!adj->is_visible() || adj->get_rel_dist_to_camera() > DRAW_DIST_TILES) continue;
			ix_sz_pair const &ixsz(crack_ibuf.lookup(crack_ibuf.get_index(dim, dir, lod_level, adj->get_lod_level(reflection_pass))));
			if (ixsz.sz == 0) continue;
			glDrawRangeElements(GL_TRIANGLES, 0, stride*stride, ixsz.sz, GL_UNSIGNED_INT, (void *)(ixsz.ix*sizeof(unsigned)));
		} // for dir
	} // for dim
	bind_vbo(0, 0); // unbind vertex buffer
	fgPopMatrix();
	if (draw_near_water) {draw_water_cap(s, 1);}
}


void tile_t::draw_water_cap(shader_t &s, bool textures_already_set) const {

	bool const dist_water(draw_distant_water());
	//if (dist_water) return; // skip water cap
	cube_t const bcube(get_mesh_bcube());
	static vector<vert_wrap_t> wverts;
	wverts.resize(0);

	// draw vertical edges that cap the water volume and will be blended between underwater black and fog colors
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			int dxy[2] = {0,0};
			dxy[dim] = (dir ? 1 : -1);
			tile_t const *const adj_tile(get_adj_tile(dxy[0], dxy[1]));
				
			if (adj_tile && !adj_tile->is_distant) {
				if (adj_tile->get_rel_dist_to_camera() < DRAW_DIST_TILES) continue;
				if (!adj_tile->has_water ()) continue; // adj tile has no water, so we can't have any uncapped water on this edge
				if (!adj_tile->is_visible()) continue; // adj tile not visible,  so we can't see this edge
			}
			float const ztop(water_plane_z - (dist_water ? 0.1 : 0.0));
			float const x1(bcube.d[0][dim ? 0 : dir]), x2(bcube.d[0][dim ? 1 : dir]), y1(bcube.d[1][dim ? dir : 0]), y2(bcube.d[1][dim ? dir : 1]);
			wverts.push_back(vert_wrap_t(point(x1, y1, mzmin)));
			wverts.push_back(vert_wrap_t(point(x2, y2, mzmin)));
			wverts.push_back(vert_wrap_t(point(x2, y2, ztop)));
			wverts.push_back(vert_wrap_t(point(x1, y1, ztop)));
		}
	}
	if (!wverts.empty()) {
		if (!textures_already_set) {bind_textures();}
		int const loc(s.get_uniform_loc("htex_scale"));
		if (loc >= 0) {s.set_uniform_float(loc, 0.0);} // disable height texture
		draw_quad_verts_as_tris(wverts);
		if (loc >= 0) {s.set_uniform_float(loc, 1.0);} // enable height texture
	}
}


bool tile_t::is_water_visible() const {
	return (!is_distant && has_water() && get_rel_dist_to_camera() < DRAW_DIST_TILES && is_visible());
}
void tile_t::draw_water(shader_t &s, float z) const {

	if (!is_water_visible()) return;
	float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2)), xv2(xv1+(x2-x1)*deltax), yv2(yv1+(y2-y1)*deltay);
	bind_texture_tu(height_tid, 2);
	draw_one_tquad(xv1, yv1, xv2, yv2, z, (use_water_plane_tess() ? GL_PATCHES : GL_TRIANGLE_FAN));
}


bool tile_t::check_sphere_collision(point &pos, float radius) const {

	if (is_distant || !contains_point(pos)) return 0;
	if (pos.z > get_tile_zmax() + radius)   return 0; // sphere is completely above the tile
	bool coll(0);

	if (!pine_trees.empty()) {
		pos  -= ptree_off.get_xlate();
		coll |= pine_trees.check_sphere_coll(pos, radius);
		pos  += ptree_off.get_xlate();
	}
	if (!decid_trees.empty()) {
		pos  -= dtree_off.get_xlate();
		coll |= decid_trees.check_sphere_coll(pos, radius);
		pos  += dtree_off.get_xlate();
	}
	if (scenery.generated) {
		pos  -= scenery_off.get_xlate();
		coll |= scenery.check_sphere_coll(pos, radius);
		pos  += scenery_off.get_xlate();
	}
	return coll;
}


int tile_t::get_tid_under_point(point const &pos) const {

	if (is_distant || !contains_point(pos)) return -1;
	int const xpos(max(0, min((int)size, (get_xpos(pos.x) - x1 - xoff + xoff2)))); // min/max not needed?
	int const ypos(max(0, min((int)size, (get_ypos(pos.y) - y1 - yoff + yoff2))));
	unsigned const ix(4*(ypos*stride + xpos));
	assert(ix < weight_data.size());
	unsigned max_weight(0), max_weight_ix(0), weight_sum(0);

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const w(weight_data[ix + i]);
		weight_sum += w;
		if (w > max_weight) {max_weight = w; max_weight_ix = i;}
	}
	if ((255 - weight_sum) > max_weight) {max_weight_ix = 4;} // 5th weight (255 - sum(w)) is max
	return lttex_dirt[max_weight_ix].id;
}


bool tile_t::line_intersect_mesh(point const &v1, point const &v2, float &t, int &xpos, int &ypos) const {

	if (is_distant) return 0; // Note: this can be made to work, but won't work as-is
	point v1c(v1), v2c(v2); // clipped verts
	if (!do_line_clip(v1c, v2c, get_mesh_bcube().d)) return 0;
	// similar to mesh_intersector::line_intersect_surface_fast()
	int const xp1(get_xpos(v1c.x) - x1 - xoff + xoff2), yp1(get_ypos(v1c.y) - y1 - yoff + yoff2);
	int const xp2(get_xpos(v2c.x) - x1 - xoff + xoff2), yp2(get_ypos(v2c.y) - y1 - yoff + yoff2);
	int const dx(xp2 - xp1), dy(yp2 - yp1), steps(max(1, max(abs(dx), abs(dy))));
	double const dz(v2c.z - v1c.z), xinc(dx/(double)steps), yinc(dy/(double)steps), zinc(dz/(double)steps);
	double x(xp1), y(yp1), z(v1c.z - 0.1*fabs(zinc)); // z offset required to avoid problems with zval at bcube.z1

	for (int k = 0; k <= steps; ++k) {
		int const ix((int)x), iy((int)y);

		if (ix >= 0 && iy >= 0 && ix <= (int)size && iy <= (int)size && zvals[iy*zvsize + ix] > z) {
			// Note: we use z instead of zvals here because zvals may be much too high if we enter this tile while the line is under the mesh
			float const cur_t(((z - 0.5*zinc) - v1.z)/(v2.z - v1.z)); // t relative to original v1, v2

			if (cur_t >= 0.0 && cur_t <= 1.0) {
				xpos = x1 + ix;
				ypos = y1 + iy;
				t    = cur_t;
				return 1;
			}
		}
		x += xinc; y += yinc; z += zinc;
	}
	return 0;
}


// *** lightning_strike_t ***

point lightning_strike_t::get_pos() const {

	assert(path.points.size() >= 2);
	//return path.points[path.points.size()-2];
	point avg_pos(all_zeros);
	for (vector<point>::const_iterator i = path.points.begin(); i != path.points.end(); ++i) {avg_pos += *i;}
	return avg_pos / path.points.size();
}

void lightning_strike_t::gen() {

	float const cloud_zmax(get_cloud_zmax()), hmax(cloud_zmax - water_plane_z);
	float const gen_radius(LITNING_DIST*get_scaled_tile_radius()), max_dxy(0.05*hmax), max_dz(0.1*hmax);
	point const camera(get_camera_pos());
	point pos((camera.x + gen_radius*signed_rand_float()), (camera.y + gen_radius*signed_rand_float()), cloud_zmax);
	path.points.push_back(pos);
		
	while (1) {
		pos.z -= max_dz*rand_float(); // move down
		for (unsigned d = 0; d < 2; ++d) {pos[d] += max_dxy*signed_rand_float();}
		float const mzval(interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 0)), zstop(max(mzval, water_plane_z));
		path.points.push_back(pos);
		if (pos.z < zstop) break; // mesh or water intersection (should the line be clipped?)
	}
	assert(path.points.size() >= 2);
	path.width = 10.0;
	path.color = LITN_C;
	time       = 1; // nonzero

	// play thunder sound
	point const lightning_pos(get_pos());
	float const delay(2.0*p2p_dist(camera, lightning_pos)/gen_radius); // max of 2s delay
	play_thunder(lightning_pos, 10.0, delay); // second to the last (above the mesh)
}

void lightning_strike_t::update() {

	int const rnum(int(fticks*LIGHTNING_FREQ));

	if (enabled()) {
		if (time > int(fticks*LITNING_TIME2)) { // check for end time
			clear();
		}
		else if (animate2) {
			time += iticks;
		}
	}
	else if (animate2 && (rnum <= 1 || (rand()%rnum) == 0)) {
		gen();
	}
}

void lightning_strike_t::draw() const {

	if (!enabled()) return;
	path.draw_lines(); // uses a custom shader with no fog
	colorRGBA const ambient(path.color*0.2);
	float const radius(0.4*get_scaled_tile_radius());
	set_colors_and_enable_light(LIGHTNING_LIGHT, ambient, path.color);
	setup_gl_light_atten(LIGHTNING_LIGHT, 0.1, 0.0, 1.0/(radius*radius));
	set_gl_light_pos(LIGHTNING_LIGHT, get_pos(), 1.0); // point light source position
	tt_lightning_enabled = 1;
}

void lightning_strike_t::end_draw() const {

	disable_light(LIGHTNING_LIGHT); // even if not currently enabled, in case it was enabled before an update
	tt_lightning_enabled = 0;
}


// *** tile_draw_t ***


tile_draw_t::tile_draw_t() : lod_renderer(USE_TREE_BILLBOARDS), terrain_zmin(0.0) {

	assert(MESH_X_SIZE == MESH_Y_SIZE && X_SCENE_SIZE == Y_SCENE_SIZE);
}


void tile_draw_t::clear() {

	clear_vbos_tids(); // needed to clear vbo, ivbo, and free list
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear();} // may not be necessary
	to_draw.clear();
	tiles.clear();
	tree_trunk_pts.clear();
}


float tile_draw_t::update(float &min_camera_dist) { // view-independent updates; returns terrain zmin

	//RESET_TIME;
	if (terrain_hmap_manager.maybe_load(mh_filename_tt, (invert_mh_image != 0))) {
		read_default_hmap_modmap();
	}
	to_draw.clear();
	terrain_zmin = FAR_DISTANCE;
	grass_tile_manager.update(); // every frame, even if not in tiled terrain mode?
	assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
	point const cpos(get_camera_pos()), camera(cpos - point((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
	int const tile_radius(int(CREATE_DIST_TILES*TILE_RADIUS) + 1);
	int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
	int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
	int const x2( tile_radius + toffx), y2( tile_radius + toffy);
	unsigned const init_tiles((unsigned)tiles.size());
	unsigned num_erased(0);
	min_camera_dist = FAR_DISTANCE;
	// Note: we may want to calculate distant low-res or larger tiles when the camera is high above the mesh

	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ) { // update tiles and free old tiles (Note: no ++i)
		if (!i->second->update_range(smap_manager)) {
			i->second->clear();
			tiles.erase(i++);
			++num_erased;
		}
		else {++i;}
	}
	for (int y = y1; y <= y2; ++y ) { // create new tiles
		for (int x = x1; x <= x2; ++x ) {
			tile_xy_pair const txy(x, y);
			if (tiles.find(txy) != tiles.end()) continue; // already exists
			tile_t tile(get_tile_size(), x, y);
			if (tile.get_rel_dist_to_camera() >= CREATE_DIST_TILES) continue; // too far away to create
			tile_t *new_tile(new tile_t(tile));
			new_tile->create_zvals(height_gen);
			tiles[txy].reset(new_tile);
		}
	}
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) { // calculate terrain_zmin
		if (i->second->get_rel_dist_to_camera() <= DRAW_DIST_TILES) {
			terrain_zmin = min(terrain_zmin, i->second->get_zmin());
			if (!camera_surf_collide) {min_camera_dist = min(min_camera_dist, i->second->get_min_dist_to_pt(cpos, 0, 0));}
		}
	}
	if (DEBUG_TILES && (tiles.size() != init_tiles || num_erased > 0)) {
		cout << "update: tiles: " << init_tiles << " to " << tiles.size() << ", erased: " << num_erased << endl;
	}
	// Note: could skip shadow computation (but not weight calc/texture upload) if (max(sun_pos.z,  last_sun.z) > zbottom) or (sun.get_norm().z > 0.9) or something like that
	static point last_sun(all_zeros), last_moon(all_zeros);
	bool const sun_change (sun_pos  != last_sun  && light_factor >= 0.4);
	bool const moon_change(moon_pos != last_moon && light_factor <= 0.6);

	// Note: we could regen trees and scenery if water was just turned on to remove underwater vegetation
	if (mesh_shadows_enabled() && (sun_change || moon_change)) { // light source change
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->clear_shadows();
		}
		last_sun  = sun_pos;
		last_moon = moon_pos;
	}
	//if ((GET_TIME_MS() - timer1) > 100) {PRINT_TIME("Tiled Terrain Update");}
	return terrain_zmin;
}

float tile_draw_t::get_actual_zmin() const {return min(zmin, terrain_zmin);}


void tile_draw_t::setup_terrain_textures(shader_t &s, unsigned start_tu_id) {

	unsigned const base_tsize(NORM_TEXELS);

	for (int i = 0; i < NTEX_DIRT; ++i) {
		int const tid(lttex_dirt[i].id);
		float const tscale(float(base_tsize)/float(get_texture_size(tid, 0))); // assumes textures are square
		float cscale(1.0);
		if (tid == GROUND_TEX) {cscale = TT_GRASS_COLOR_SCALE;} // darker grass
		if (tid == ROCK_TEX  ) {cscale = 0.5;} // darker rock
		unsigned const tu_id(start_tu_id + i);
		select_multitex(tid, tu_id);
		std::ostringstream oss1, oss2, oss3;
		oss1 << "tex" << tu_id;
		oss2 << "ts"  << tu_id;
		oss3 << "cs"  << tu_id;
		s.add_uniform_int(  oss1.str().c_str(), tu_id);
		s.add_uniform_float(oss2.str().c_str(), tscale);
		s.add_uniform_float(oss3.str().c_str(), cscale);
	}
	s.add_uniform_color("snow_cscale", colorRGB(1.0, 1.0, 1.2)); // increased blue
}


void setup_tt_fog_pre(shader_t &s) {

	s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
	if (nonunif_fog_enabled()) {s.set_prefix("#define USE_NONUNIFORM_FOG", 1);} // FS
	if (volume_lighting && is_cloudy && light_factor >= 0.6) {s.set_prefix("#define GOD_RAYS", 1);} // FS (sun out and cloudy)
}

void setup_tt_fog_post(shader_t &s) {

	s.setup_fog_scale();
	s.add_uniform_float("fog_bot", get_tt_fog_top());
	s.add_uniform_float("fog_top", get_tt_fog_bot());
	s.add_uniform_float("eye_z",   get_camera_pos().z);
}

void tile_draw_t::shared_shader_lighting_setup(shader_t &s, unsigned lighting_shader) {

	s.setup_enabled_lights(3, (1 << lighting_shader)); // sun, moon, and lightning
	if (!underwater && clouds_enabled()) {s.set_prefix("#define FOG_FADE_TO_TRANSPARENT", 1);} // FS - fade distant hills into background clouds
	setup_tt_fog_pre(s);
}

void tile_draw_t::lighting_with_cloud_shadows_setup(shader_t &s, unsigned lighting_shader, bool cloud_shadows) {

	shared_shader_lighting_setup(s, lighting_shader);
	s.set_prefix("#define NUM_OCTAVES 4", lighting_shader); // for clouds
	s.set_bool_prefix("apply_cloud_shadows", cloud_shadows, lighting_shader); // FS
}

void setup_cloud_plane_uniforms(shader_t &s, float cloud_cover_factor=0.535, bool match_cloud_layer=0) {

	float cloud_zmax;
	if (match_cloud_layer) {cloud_zmax = get_cloud_zmax();} // follows the camera zval - matches the drawn cloud layer but moves clouds on the terrain
	else {cloud_zmax = 0.5*(zmin + zmax) + max(zmax, CLOUD_CEILING);} // fixed z value - independent of camera z so stays in place, but disagrees with drawn clouds
	set_cloud_uniforms(s, 9);
	s.add_uniform_float("cloud_scale",   (is_cloudy ? 1.0 : (cloud_cover + cloud_cover_factor*(1.0 - cloud_cover))));
	s.add_uniform_float("cloud_alpha",   (is_cloudy ? 0.8 : 0.75)*atmosphere);
	s.add_uniform_float("cloud_plane_z", cloud_zmax);
}

void set_tile_xy_vals(shader_t &s) {

	float const inv_scale(1.0/(get_tile_size() + 1.0));
	s.add_uniform_float("x1", -0.5*DX_VAL);
	s.add_uniform_float("y1", -0.5*DY_VAL);
	s.add_uniform_float("dx_inv", inv_scale*DX_VAL_INV);
	s.add_uniform_float("dy_inv", inv_scale*DY_VAL_INV);
}


// uses texture units 0-11 (12 if using hmap texture, 13-14 if using shadow maps)
// Note: could be static, except uses get_actual_zmin()
void tile_draw_t::setup_mesh_draw_shaders(shader_t &s, bool reflection_pass, bool enable_shadow_map) const {

	bool const has_water(is_water_enabled() && !reflection_pass);
	lighting_with_cloud_shadows_setup(s, 1, (cloud_shadows_enabled() && !reflection_pass));
	bool const water_caustics(has_water && !(display_mode & 0x80) && (display_mode & 0x100) && water_params.alpha < 1.5);
	bool const use_normal_map(!reflection_pass && (display_mode & 0x08) != 0); // enabled by default
	bool const rain_mode(is_precip_enabled() && temperature > W_FREEZE_POINT);
	if (has_water      ) {s.set_prefix("#define HAS_WATER",       1);} // FS
	if (water_caustics ) {s.set_prefix("#define WATER_CAUSTICS",  1);} // FS
	if (use_normal_map ) {s.set_prefix("#define USE_NORMAL_MAP",  1);} // FS
	if (reflection_pass) {s.set_prefix("#define REFLECTION_MODE", 1);} // FS
	s.set_bool_prefix("use_shadow_map", enable_shadow_map, 1); // FS
	s.set_vert_shader("tiled_mesh");
	s.set_frag_shader("water_fog.part*+linear_fog.part+perlin_clouds.part*+ads_lighting.part*+shadow_map.part*+detail_normal_map.part+triplanar_texture.part+tiled_mesh");
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("weights_tex", 0);
	s.add_uniform_int("detail_tex",  1);
	s.add_uniform_int("shadow_normal_tex", 7);
	s.add_uniform_float("normal_z_scale", (reflection_pass ? -1.0 : 1.0));
	s.add_uniform_float("spec_offset", (rain_mode ? 0.5 : 0.0)); // increase specular during rain
	if (enable_shadow_map) {s.add_uniform_float("smap_atten_cutoff", get_smap_atten_val());}
	set_noise_tex(s, 8);
	setup_cloud_plane_uniforms(s);
	setup_terrain_textures(s, 2);
	set_tile_xy_vals(s);
	s.add_uniform_int("height_tex", 12);
	if (use_normal_map) {setup_detail_normal_map(s, 2.0);}
	
	if (enable_shadow_map) {
		s.add_uniform_float("z_bias", DEF_Z_BIAS);
		s.add_uniform_float("pcf_offset", shadow_map_pcf_offset);
	}
	if (has_water) {
		set_water_plane_uniforms(s);
		s.add_uniform_float("water_atten",    WATER_COL_ATTEN*mesh_scale);
		s.add_uniform_color("uw_atten_max",   uw_atten_max);
		s.add_uniform_color("uw_atten_scale", uw_atten_scale);

		if (water_caustics) {
			select_multitex(WATER_CAUSTIC_TEX, 10);
			s.add_uniform_int("caustic_tex", 10);
			s.add_uniform_float("caustics_weight", (1.5 - water_params.alpha));
		}
		setup_water_plane_texgen(8.0, 2.5, s, 2); // increase texture scale and change AR since the caustics texture is sparser than the water texture
	}
	else { // or just disable water fog calculation in the vertex shader (water_fog.part)?
		s.add_uniform_float("water_plane_z", (reflection_pass ? water_plane_z : get_actual_zmin())); // used for fog calculation/clipping
	}
	set_landscape_texgen(1.0, -MESH_X_SIZE/2, -MESH_Y_SIZE/2, MESH_X_SIZE, MESH_Y_SIZE, s, 1);
	s.add_uniform_float("triplanar_texture_scale", 1.0/(X_SCENE_SIZE + Y_SCENE_SIZE));
}


bool tile_draw_t::can_have_reflection_recur(tile_t const *const tile, point const corners[3], tile_set_t &tile_set, unsigned dim_ix) {

	point const camera(get_camera_pos());
	cube_t bcube(tile->get_bcube());
	bcube.d[2][0] = min(bcube.d[2][0], zmin); // make sure the line isn't clipped in z
	bcube.d[2][1] = max(bcube.d[2][1], zmax);
	if (dim_ix < 2 && !check_line_clip(camera, corners[dim_ix], bcube.d)) return 0; // not within the shadow of the original tile

	if (!tile->was_last_occluded() && tile->has_water()) {
		cube_t water_bcube(tile->get_water_bcube());

		if (camera_pdu.cube_visible(water_bcube)) {
			point const closest_pt(water_bcube.closest_pt(corners[2]));
			float const z_over_xy((camera.z - water_plane_z)/p2p_dist_xy(camera, closest_pt)); // slope
			if (water_plane_z + z_over_xy*p2p_dist_xy(corners[2], closest_pt) < corners[2].z) return 1;
		}
	}
	if (!tile_set.insert(tile->get_tile_xy_pair()).second) return 0; // already seen
		
	for (unsigned d = 0; d < 2; ++d) {
		int delta[2] = {0,0};
		if      (camera[d] < bcube.d[d][0]) {delta[d] = -1;}
		else if (camera[d] > bcube.d[d][1]) {delta[d] =  1;}
		else                                {continue;}
		tile_t const *const adj(get_tile_from_xy(tile->get_tile_xy_pair(delta[0], delta[1])));
		if (!adj || !adj->is_visible()) continue;
		if (can_have_reflection_recur(adj, corners, tile_set, d)) return 1;
	}
	return 0;
}


bool tile_draw_t::can_have_reflection(tile_t const *const tile, tile_set_t &tile_set) {

	if (tile->all_water()) return 0;
	if (tile->has_water()) return 1;
	// return 1 if the camera's tile contains water?
	point const camera(get_camera_pos());
	cube_t bcube(tile->get_bcube());
	point const center(bcube.get_cube_center());
	point corners[3]; // {x, y, closest}

	for (unsigned d = 0; d < 2; ++d) {
		corners[ d][d] = bcube.d[d][center[d] < camera[d]]; // 0,0
		corners[!d][d] = bcube.d[d][center[d] > camera[d]]; // 1,1
		corners[ d][2] = bcube.d[2][0];
	}
	corners[2]   = bcube.closest_pt(camera);
	corners[2].z = tile->get_zmax();
	return can_have_reflection_recur(tile, corners, tile_set, 2);
}


void tile_draw_t::pre_draw() { // view-dependent updates/GPU uploads

	vector<tile_t *> to_update, to_gen_trees;
	assert((vbo == 0) == (ivbo == 0)); // either neither or both are valid
	
	if (vbo == 0) { // build mesh vbo/ivbo
		unsigned const tile_size(get_tile_size()), stride(tile_size+1);
		vector<point> data(stride*stride);
		vector<unsigned> indices;

		for (unsigned y = 0; y <= tile_size; ++y) {
			for (unsigned x = 0; x <= tile_size; ++x) {
				data[y*stride + x].assign(x*DX_VAL, y*DY_VAL, 0.0); // z=0.0
			}
		}
		for (unsigned i = 0, step = 1; i < NUM_LODS; ++i, step <<= 1) {
			ivbo_ixs[i] = indices.size();
			unsigned const isz_ceil((tile_size + step - 1)/step);
			if (isz_ceil < 2) continue; // too small, don't create LOD for this level (will never be used during rendering)

			for (unsigned y = 0, ny = 0; ny < isz_ceil; y += step, ++ny) {
				for (unsigned x = 0, nx = 0; nx <= isz_ceil; x += step, ++nx) { // 2 extra to start the strip
					indices.push_back(y*stride + x);
					indices.push_back(min(y+step, tile_size)*stride + x);
				}
				indices.push_back(PRIMITIVE_RESTART_IX); // restart the strip
			}
		} // for i
		ivbo_ixs[NUM_LODS] = indices.size();
		crack_ibuf.gen_offsets(indices, tile_size);
		create_and_upload(data, indices, 0, 1); // unbind at end
	}
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		tile_t *const tile(i->second.get());
		assert(tile);
		if (tile->get_rel_dist_to_camera() > DRAW_DIST_TILES) continue; // too far to draw
		if (!tile->is_visible()) continue; // Note: using current camera view frustum

		if (tile->can_have_trees()) { // no trees in water or distant tiles
			if (pine_trees_enabled() && !tile->pine_trees_generated()) {to_gen_trees.push_back(tile);}
			if (decid_trees_enabled()) {tile->gen_decid_trees_if_needed();}
		}
		to_update.push_back(tile);
	}
	if (enable_instanced_pine_trees() && !to_gen_trees.empty()) {
		create_pine_tree_instances();
	}
	//RESET_TIME;
	// don't use parallel tree gen for a single tile, or when GPU heightmaps are enabled
	#pragma omp parallel for schedule(dynamic,1) if (mesh_gen_mode != 3 && to_gen_trees.size() > 1)
	for (int i = 0; i < (int)to_gen_trees.size(); ++i) {to_gen_trees[i]->init_pine_tree_draw();}
	//PRINT_TIME("Gen Trees2");
	
	for (vector<tile_t *>::iterator i = to_update.begin(); i != to_update.end(); ++i) {
		(*i)->pre_draw(height_gen);

		if ((*i)->can_have_trees()) {
			(*i)->update_pine_tree_state(1);
			(*i)->update_decid_trees();
		}
		(*i)->update_scenery();
	}
	for (vector<tile_t *>::iterator i = to_update.begin(); i != to_update.end(); ++i) { // after everything has been setup
		(*i)->setup_shadow_maps(smap_manager);
	}
}


void tile_draw_t::occluder_pts_t::calc_cube_top_points(cube_t const &bcube) { // copied from cube_t::get_points()

	unsigned i[3] = {0,0,1};

	for (i[0] = 0; i[0] < 2; ++i[0]) {
		for (i[1] = 0; i[1] < 2; ++i[1]) {
			UNROLL_3X(cube_pts[(i[0]<<1)+i[1]][i_] = bcube.d[i_][i[i_]];)
		}
	}
}


unsigned in_mb(unsigned long long v) {return v/1024/1024;}


void tile_draw_t::draw(bool reflection_pass) {

	//cout << "zmin: " << zmin << ", zmax: " << zmax << ", wpz: " << water_plane_z << endl; // TESTING
	//RESET_TIME;
	unsigned num_trees(0);
	unsigned long long mem(0), tree_mem(0);
	to_draw.clear();
	occluded_tiles.clear();

	if (DEBUG_TILES && vbo) {
		unsigned const tile_size(get_tile_size());
		mem += 2*tile_size*(tile_size+1)*sizeof(point);
		for (unsigned i = 0; i < NUM_LODS; ++i) {mem += 4*(tile_size>>i)*(tile_size>>i)*sizeof(unsigned);} // approximate
	}

	// determine potential occluders
	vector<tile_t *> occluders;
	vector<cube_t> test_cubes;
	point const camera(get_camera_pos());

	if ((display_mode & 0x08) && (display_mode & 0x01)) { // check occlusion when occlusion culling and mesh are enabled
		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
			tile_t *const tile(i->second.get());
			if (tile->get_rel_dist_to_camera() > OCCLUDER_DIST || !tile->is_visible()) continue;
			occluders.push_back(tile);
		}
	}
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		tile_t *const tile(i->second.get());
		if (DEBUG_TILES) {mem      += tile->get_gpu_mem ();}
		if (DEBUG_TILES) {tree_mem += tile->get_tree_mem();}
		float const dist(tile->get_rel_dist_to_camera());
		if (dist > DRAW_DIST_TILES || !tile->is_visible()) continue;
		if (tile->was_last_occluded()) continue; // occluded in the shadow pass
		tile_set_t tile_set;
		if (reflection_pass && !can_have_reflection(tile, tile_set)) continue;

		if (!occluders.empty() && !tile->was_last_unoccluded()) {
			occluder_pts_t tile_os, sub_tile_os;
			tile_os.calc_cube_top_points(tile->get_bcube());
			bool tile_occluded(1);
			test_cubes.resize(0);

			for (vector<tile_t *>::const_iterator j = occluders.begin(); j != occluders.end(); ++j) {
				if (*j == tile) continue; // no self-occlusion
				cube_t occluder_bcube((*j)->get_mesh_bcube());
				occluder_bcube.d[2][0] = zmin; // not required?
				bool intersected(0); // will always be true for the camera's current tile

				for (unsigned d = 0; d < 4 && !intersected; ++d) {
					intersected |= check_line_clip(tile_os.cube_pts[d], camera, occluder_bcube.d);
				}
				if (!intersected) continue; // skip

				for (unsigned s = 0; s < 16; ++s) {
					test_cubes.push_back((*j)->get_mesh_sub_bcube((s>>2), (s&3)));
					test_cubes.back().d[2][1] = test_cubes.back().d[2][0]; // cube below the bcube
					test_cubes.back().d[2][0] = zmin;
				}
			} // for j
			for (unsigned t = 0; t < 16; ++t) {
				sub_tile_os.calc_cube_top_points(tile->get_mesh_sub_bcube((t>>2), (t&3)));
				bool sub_tile_occluded(0);

				for (vector<cube_t>::const_iterator s = test_cubes.begin(); s != test_cubes.end() && !sub_tile_occluded; ++s) {
					sub_tile_occluded = 1;

					for (unsigned d = 0; d < 4 && sub_tile_occluded; ++d) {
						sub_tile_occluded &= check_line_clip(sub_tile_os.cube_pts[d], camera, s->d);
					}
				}
				if (!sub_tile_occluded) {tile_occluded = 0; break;}
			} // for t
			tile->set_last_occluded(tile_occluded);
			if (tile_occluded) {occluded_tiles.push_back(tile); continue;}
		} // check_occlusion
		to_draw.push_back(make_pair(dist, tile));
		num_trees += tile->num_pine_trees() + tile->num_decid_trees();
	} // for i
	sort(to_draw.begin(), to_draw.end()); // sort front to back to improve draw time through depth culling

	if (display_mode & 0x01) { // draw visible tiles
		if (shadow_map_enabled()) {draw_tiles(reflection_pass, 1);} // shadow map pass
		draw_tiles(reflection_pass, 0); // non-shadow map pass
	}
	
	if (DEBUG_TILES) {
		unsigned const dtree_mem(tree_data_manager.get_gpu_mem()), ptree_mem(get_pine_tree_inst_gpu_mem()), grass_mem(grass_tile_manager.get_gpu_mem());
		cout << "tiles drawn: " << to_draw.size() << " of " << tiles.size()
			 << ", trees drawn: " << num_trees << ", gpu mem: " << in_mb(mem + tree_mem + dtree_mem + ptree_mem + grass_mem)
			 << ", tree mem: " << in_mb(tree_mem) << ", decid tree mem: " << in_mb(dtree_mem) << ", grass mem: " << in_mb(grass_mem) << endl;
	}
	if (pine_trees_enabled ()) {draw_pine_trees (reflection_pass);}
	if (decid_trees_enabled()) {draw_decid_trees(reflection_pass);}
	if (scenery_enabled    ()) {draw_scenery    (reflection_pass);}
	if (is_grass_enabled   ()) {draw_grass      (reflection_pass);}
	lightning_strike.end_draw(); // in case it was enabled
	//if ((GET_TIME_MS() - timer1) > 100) {PRINT_TIME("Draw Tiled Terrain");}
}


void tile_draw_t::draw_tiles(bool reflection_pass, bool enable_shadow_map) const {

	shader_t s;
	setup_mesh_draw_shaders(s, reflection_pass, enable_shadow_map);
	s.add_uniform_float("spec_scale", 1.0);
	s.clear_specular(); // in case we failed to clear it somewhere ahead
	s.enable_vnct_atribs(1, 0, 0, 0);
	enable_blend(); // for fog transparency
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(PRIMITIVE_RESTART_IX);

	for (unsigned i = 0; i < to_draw.size(); ++i) {
		if (to_draw[i].second->using_shadow_maps() != enable_shadow_map) continue; // draw in another pass
		to_draw[i].second->draw(s, *this, ivbo_ixs, crack_ibuf, reflection_pass);
	}
	if (!enable_shadow_map && !reflection_pass && is_water_enabled()) { // draw all water caps (required, even for occluded tiles)
		for (auto i = occluded_tiles.begin(); i != occluded_tiles.end(); ++i) {(*i)->draw_water_cap(s, 0);}
	}
	bind_vbo(0, 1); // unbind index buffer
	glDisable(GL_PRIMITIVE_RESTART);
	disable_blend();

	if ((display_mode & 0x01) && draw_distant_water()) {
		int const loc(s.get_uniform_loc("htex_scale"));
		if (loc >= 0) {s.set_uniform_float(loc, 0.0);} // disable height texture
		draw_distant_mesh_bottom(terrain_zmin); // Note: textures from last drawn tile are bound, but don't really affect the results
		if (loc >= 0) {s.set_uniform_float(loc, 1.0);} // enable height texture
	}
	s.end_shader();
}


void tile_draw_t::draw_shadow_pass(point const &lpos, tile_t *tile) {

	//RESET_TIME;
	bool const orig_fog_enabled(fog_enabled);
	fog_enabled = 0; // optimization?
	to_draw.clear();

	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		if (!i->second->is_visible()) continue;
		i->second->update_pine_tree_state(1);
		//i->second->update_decid_trees(); // not legal
		to_draw.push_back(make_pair(0.0, i->second.get())); // distance is unused so set to 0.0
	}
	if (pine_trees_enabled ()) {draw_pine_trees (0, 1);}
	if (decid_trees_enabled()) {draw_decid_trees(0, 1);}
	if (scenery_enabled    ()) {draw_scenery    (0);}
	render_models(1, model3d_offset.get_xlate()); // VFC should work here (somewhat?) for models
	fog_enabled = orig_fog_enabled;
	//PRINT_TIME("Draw Shadow Pass");
}


void tile_draw_t::draw_water(shader_t &s, float zval) const {
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->draw_water(s, zval);}
}


void tile_draw_t::set_noise_tex(shader_t &s, unsigned tu_id) {

	select_multitex(DITHER_NOISE_TEX, tu_id, 1);
	s.add_uniform_int("noise_tex", tu_id);
}

void tile_draw_t::set_tree_dither_noise_tex(shader_t &s, unsigned tu_id) {

	set_noise_tex(s, tu_id);
	s.add_uniform_float("noise_tex_size", get_texture_size(DITHER_NOISE_TEX, 0));
}

void tile_draw_t::set_pine_tree_shader(shader_t &s, string const &vs) {

	shared_shader_lighting_setup(s, 0);
	s.set_vert_shader("ads_lighting.part*+texture_gen.part+" + vs);
	s.set_frag_shader("linear_fog.part+noise_dither.part+pine_tree");
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("branch_tex", 0);
	s.add_uniform_float("min_alpha", 0.5);
	set_tree_dither_noise_tex(s, 1); // TU=1
	check_gl_error(302);
}


void tile_draw_t::draw_pine_tree_bl(shader_t &s, bool branches, bool near_leaves, bool far_leaves, bool reflection_pass, int xlate_loc) {

	for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
		to_draw[i].second->draw_pine_trees(s, tree_trunk_pts, branches, near_leaves, far_leaves, reflection_pass, xlate_loc);
	}
}


void tile_draw_t::draw_pine_trees(bool reflection_pass, bool shadow_pass) {

	// far leaves
	if (!shadow_pass) {
		enable_blend(); // for fog transparency
		shader_t s;
		set_pine_tree_shader(s, "pine_tree_billboard_auto_orient");
		s.add_uniform_float("radius_scale", calc_tree_size());
		s.add_uniform_float("ambient_scale", 1.5);
		s.set_specular(0.2, 8.0);
		draw_pine_tree_bl(s, 0, 0, 1, reflection_pass, s.get_uniform_loc("xlate"));
		s.clear_specular();
		s.end_shader();
		disable_blend();
	}
	// near leaves
	shader_t s;
	if (shadow_pass) {s.set_prefix("#define NO_NOISE", 1);} // FS
	if (enable_instanced_pine_trees()) {s.set_prefix("#define ENABLE_INSTANCING", 0);} // VS
	set_pine_tree_shader(s, "pine_tree");
	int xlate_loc(-1);
	
	if (enable_instanced_pine_trees()) {
		s.add_uniform_float("vertex_scale", calc_tree_size()); // default is 0.8
		xlate_loc = s.get_attrib_loc("xlate");
		glEnableVertexAttribArray(xlate_loc);
		glVertexAttribDivisor(xlate_loc, 1);
	}
	s.set_specular(0.2, 8.0);
	draw_pine_tree_bl(s, 0, 1, 0, reflection_pass, xlate_loc);
	assert(tree_trunk_pts.empty());
	s.clear_specular();
	
	if (xlate_loc >= 0) {
		glVertexAttribDivisor(xlate_loc, 0);
		glDisableVertexAttribArray(xlate_loc);
	}
	s.end_shader();

	// nearby trunks
	setup_tt_fog_pre(s);
	setup_smoke_shaders(s, 0.0, 0, 0, 0, !shadow_pass, 0, 0, 0, 0, 0, 0, 1);
	setup_tt_fog_post(s);
	s.add_uniform_color("const_indir_color", colorRGB(0,0,0)); // don't want indir lighting for tree trunks
	s.add_uniform_float("tex_scale_t", 5.0);
	draw_pine_tree_bl(s, 1, 0, 0, reflection_pass); // branches
	s.add_uniform_float("tex_scale_t", 1.0);
	s.end_shader();

	// distant trunks
	if (!shadow_pass && !tree_trunk_pts.empty()) { // color/texture already set above
		enable_blend(); // for fog transparency
		set_pine_tree_shader(s, "xy_billboard");
		assert(!(tree_trunk_pts.size() & 1));
		select_texture(WHITE_TEX);
		s.set_cur_color(get_tree_trunk_color(T_PINE, 1));
		draw_verts(tree_trunk_pts, GL_LINES);
		s.end_shader();
		disable_blend();
	}
	tree_trunk_pts.resize(0);
}


void tile_draw_t::draw_decid_tree_bl(shader_t &s, tree_lod_render_t &lod_renderer, bool branches, bool leaves, bool reflection_pass) {

	for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
		to_draw[i].second->draw_decid_trees(s, lod_renderer, branches, leaves, reflection_pass);
	}
}


void tile_draw_t::billboard_tree_shader_setup(shader_t &s) {

	shared_shader_lighting_setup(s, 1);
	s.begin_shader();
	setup_tt_fog_post(s);
#ifdef USE_TREE_BB_TEX_ATLAS
	s.add_uniform_vector2d("normal_tc_off",   vector2d(0.5, 0.0));
	s.add_uniform_vector2d("normal_tc_scale", vector2d(0.5, 1.0));
	s.add_uniform_int("normal_map", 0);
#else
	s.add_uniform_vector2d("normal_tc_off",   vector2d(0.0, 0.0));
	s.add_uniform_vector2d("normal_tc_scale", vector2d(1.0, 1.0));
	s.add_uniform_int("normal_map", 1);
#endif
	s.add_uniform_int("color_map",   0);
	s.add_uniform_int("tc_start_ix", 0);
	set_tree_dither_noise_tex(s, 2); // TU=2
}


void tile_draw_t::draw_decid_trees(bool reflection_pass, bool shadow_pass) {

	float const cscale(0.8*(cloud_shadows_enabled() ? 0.75 : 1.0));
	bool const enable_billboards(USE_TREE_BILLBOARDS && !shadow_pass);
	lod_renderer.resize_zero();
	lod_renderer.set_enabled(enable_billboards); // need full detail rendering in shadow pass, since billboards project poor shadows

	{ // draw leaves
		shader_t ls;
		tree_cont_t::pre_leaf_draw(ls, enable_billboards, shadow_pass);
		if (enable_billboards) {lod_renderer.leaf_opacity_loc = ls.get_uniform_loc("opacity");}
		set_tree_dither_noise_tex(ls, 1); // TU=1
		ls.add_uniform_color("color_scale", colorRGBA(cscale, cscale, cscale, 1.0));
		draw_decid_tree_bl(ls, lod_renderer, 0, 1, reflection_pass);
		ls.add_uniform_color("color_scale", WHITE);
		tree_cont_t::post_leaf_draw(ls);
	}
	{ // draw branches
		shader_t bs;
		bs.setup_enabled_lights(3, 2); // FS; sun, moon, and lightning
		setup_tt_fog_pre(bs);
		bs.set_vert_shader("per_pixel_lighting");
		bs.set_frag_shader("linear_fog.part+ads_lighting.part*+noise_dither.part+tiled_tree_branches");
		bs.begin_shader();
		if (enable_billboards) {lod_renderer.branch_opacity_loc = bs.get_uniform_loc("opacity");}
		setup_tt_fog_post(bs);
		set_tree_dither_noise_tex(bs, 1); // TU=1
		bs.add_uniform_int("tex0", 0);
		bs.add_uniform_color("const_indir_color", colorRGB(0,0,0)); // don't want indir lighting for tree trunks
		draw_decid_tree_bl(bs, lod_renderer, 1, 0, reflection_pass);
		bs.add_uniform_vector3d("world_space_offset", zero_vector); // reset
		bs.end_shader();
	}
	lod_renderer.finalize();
	enable_blend(); // for fog transparency

	if (lod_renderer.has_leaves()) { // draw leaf billboards
		shader_t lrs;
		lrs.set_vert_shader("texture_gen.part+tree_leaves_billboard");
		lrs.set_frag_shader("linear_fog.part+ads_lighting.part*+leaf_lighting_comp.part*+noise_dither.part+tree_leaves_billboard");
		billboard_tree_shader_setup(lrs);
		lrs.add_uniform_color("color_scale", colorRGBA(cscale, cscale, cscale, 1.0));
		lrs.set_specular(0.1, 10.0);
		lod_renderer.render_billboards(0);
		lrs.clear_specular();
		lrs.end_shader();
	}
	if (lod_renderer.has_branches()) { // draw branch billboards
		shader_t brs;
		brs.set_vert_shader("texture_gen.part+tree_branches_billboard");
		brs.set_frag_shader("linear_fog.part+ads_lighting.part*+noise_dither.part+tree_branches_billboard");
		billboard_tree_shader_setup(brs); // cscale=1.0 ?
		brs.add_uniform_vector3d("ref_dir", plus_y);
		brs.add_uniform_vector3d("camera_pos", get_camera_pos());
		lod_renderer.render_billboards(1);
		brs.end_shader();
	}
	disable_blend();
	if (!reflection_pass) {leaf_color_changed = 0;} // Note: only trees in visible tiles will be updated
}


void tile_draw_t::draw_scenery(bool reflection_pass) {

	shader_t s;
	shared_shader_lighting_setup(s, 0);
	s.set_vert_shader("ads_lighting.part*+two_lights_texture");
	s.set_frag_shader("linear_fog.part+textured_with_fog");
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("tex0", 0);
		
	for (unsigned i = 0; i < to_draw.size(); ++i) {
		to_draw[i].second->draw_scenery(s, 1, 0, reflection_pass); // opaque
	}
	tree_scenery_pld.draw_and_clear();
	s.end_shader();
	set_leaf_shader(s, 0.9, 1, 0, 0, 1);

	for (unsigned i = 0; i < to_draw.size(); ++i) {
		to_draw[i].second->draw_scenery(s, 0, 1, reflection_pass); // leaves
	}
	s.end_shader();
}


void tile_draw_t::setup_grass_flower_shader(shader_t &s, bool enable_wind, bool use_smap, float dist_const_mult) {

	s.begin_shader();
	if (enable_wind) {setup_wind_for_shader(s, 1);}
	if (use_smap   ) {s.add_uniform_float("smap_atten_cutoff", get_smap_atten_val());}
	s.add_uniform_int("tex0", 0);
	s.add_uniform_int("height_tex", 2);
	s.add_uniform_int("shadow_normal_tex", 4);
	s.add_uniform_float("dist_const", dist_const_mult*get_grass_thresh());
	s.add_uniform_float("dist_slope", GRASS_DIST_SLOPE);
	setup_tt_fog_post(s);
	setup_cloud_plane_uniforms(s);
	set_tile_xy_vals(s);
}


// tu's used: 0: grass, 1: wind noise, 2: heightmap, 3: grass weight, 4: shadow map, 5: noise, 9: cloud noise
void tile_draw_t::draw_grass(bool reflection_pass) {

	if (reflection_pass) return; // no grass reflection (yet)
	bool const use_cloud_shadows(GRASS_CLOUD_SHADOWS && cloud_shadows_enabled());
	vector<vector<vector2d> > insts[NUM_GRASS_LODS];

	for (unsigned wpass = 0; wpass < 2; ++wpass) { // wind, no wind
		for (unsigned spass = 0; spass < 2; ++spass) { // shadow maps, no shadow maps
			if (spass == 0 && !shadow_map_enabled()) continue;
			shader_t s;
			bool const enable_wind((display_mode & 0x0100) && wpass == 0);
			lighting_with_cloud_shadows_setup(s, 0, use_cloud_shadows);
			if (wpass == 1) {s.set_prefix("#define DEC_HEIGHT_WHEN_FAR", 0);} // VS
			s.set_bool_prefix("enable_grass_wind", enable_wind, 0); // VS
			s.set_bool_prefix("use_shadow_map", (spass == 0), 0); // VS
			s.set_vert_shader("ads_lighting.part*+perlin_clouds.part*+shadow_map.part*+tiled_shadow_map.part*+wind.part*+grass_texture.part+grass_tiled");
			s.set_frag_shader("linear_fog.part+grass_tiled");
			//s.set_geom_shader("grass_tiled"); // triangle => triangle - too slow
			setup_grass_flower_shader(s, enable_wind, (spass == 0), 1.0);
			s.add_uniform_int("weight_tex", 3);
			set_noise_tex(s, 5);
			s.add_uniform_float("height", grass_length);
			s.set_specular(0.1, 20.0);
			grass_tile_manager.begin_draw();

			int const lt_loc(s.get_attrib_loc("local_translate"));
			glEnableVertexAttribArray(lt_loc);
			glVertexAttribDivisor(lt_loc, 1);

			for (unsigned i = 0; i < to_draw.size(); ++i) {
				if (to_draw[i].second->using_shadow_maps() != (spass == 0)) continue;
				if ((to_draw[i].second->get_dist_to_camera_in_tiles(0) > 0.5) != wpass) continue; // xyz dist
				to_draw[i].second->draw_grass(s, insts, use_cloud_shadows, lt_loc);
			}
			glVertexAttribDivisor(lt_loc, 0);
			glDisableVertexAttribArray(lt_loc);
			grass_tile_manager.end_draw();
			s.end_shader();
		} // for spass
	} // for wpass

	// draw flowers
	for (unsigned spass = 0; spass < 2; ++spass) { // shadow maps, no shadow maps
		if (flower_density == 0.0 || (display_mode & 0x10)) continue; // no flowers
		if (spass == 0 && !shadow_map_enabled()) continue;
		shader_t s;
		lighting_with_cloud_shadows_setup(s, 1, use_cloud_shadows);
		s.set_bool_prefix("use_shadow_map", (spass == 0), 1); // FS
		s.set_vert_shader("texture_gen.part+wind.part*+flowers_tiled");
		//s.set_vert_shader("texture_gen.part+wind.part*+flowers_tiled_gs");
		//s.set_geom_shader("flower_from_pt"); // point => 1 quad
		s.set_frag_shader("linear_fog.part+ads_lighting.part*+perlin_clouds.part*+shadow_map.part*+tiled_shadow_map.part*+flowers_tiled");
		setup_grass_flower_shader(s, 1, (spass == 0), FLOWER_REL_DIST);
		flower_manager_t::setup_flower_shader_post(s);
		enable_blend(); // required for distant flower transition to alpha=0

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			if (to_draw[i].second->using_shadow_maps() != (spass == 0)) continue;
			to_draw[i].second->draw_flowers(s, use_cloud_shadows);
		}
		disable_blend();
		s.end_shader();
	}
}


void tile_draw_t::update_lightning(bool reflection_pass) {

	if (!reflection_pass) {lightning_strike.update();}
	lightning_strike.draw();
}


void tile_draw_t::clear_vbos_tids() {

	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear_vbo_tid(nullptr);}
	smap_manager.clear_context();
	height_gen.clear_context();
	clear_vbos();
}

void tile_draw_t::clear_flowers() {
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear_flowers();}
}


tile_t *tile_draw_t::get_tile_from_xy(tile_xy_pair const &tp) const {
	tile_map::const_iterator it(tiles.find(tp));
	return ((it != tiles.end()) ? it->second.get() : nullptr);
}
tile_t *tile_draw_t::get_tile_containing_point(point const &pos) const {
	tile_xy_pair const tp(round_fp(0.5*(pos.x - (xoff - xoff2)*DX_VAL)/X_SCENE_SIZE), round_fp(0.5*(pos.y - (yoff - yoff2)*DY_VAL)/Y_SCENE_SIZE));
	return get_tile_from_xy(tp);
}

bool tile_draw_t::check_sphere_collision(point &pos, float radius) const { // Note: pos is modified
	tile_t const *const tile(get_tile_containing_point(pos)); // can return null for camera during tile generation frames
	return (tile ? tile->check_sphere_collision(pos, radius) : 0);
}
bool tile_draw_t::check_player_collision() const {
	point camera(get_camera_pos());
	if (!check_sphere_collision(camera, CAMERA_RADIUS)) return 0;
	surface_pos = camera; // write modified camera pos back to the scene state
	return 1;
}

int tile_draw_t::get_tid_under_point(point const &pos) const {
	tile_t const *const tile(get_tile_containing_point(pos));
	return (tile ? tile->get_tid_under_point(pos) : -1);
}


bool tile_draw_t::line_intersect_mesh(point const &v1, point const &v2, float &t, tile_t *&intersected_tile, int &xpos, int &ypos) const {

	t = 2.0; // > 1.0
	intersected_tile = nullptr;

	// Note: could start with get_tile_containing_point(v1) and walk the line to v2 using DDA, etc.
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		float tn(1.0);
		int new_xpos(0), new_ypos(0);

		// Note: could make this faster by passing tmin, tmax into line_intersect_mesh() here and using for early termination,
		// but this code is plenty fast enough to do a single query each frame as it is
		if (i->second->line_intersect_mesh(v1, v2, tn, new_xpos, new_ypos) && tn < t) {
			t = tn; xpos = new_xpos; ypos = new_ypos;
			intersected_tile = i->second.get(); // constness?
		}
	}
	if (intersected_tile != nullptr) {
		assert(t >= 0 && t <= 1.0); // okay even with fp error?
		return 1;
	}
	return 0;
}


tile_draw_t terrain_tile_draw;


void update_tiled_grass_length_width(float lscale, float wscale) {
	grass_tile_manager.scale_grass(lscale, wscale);
	terrain_tile_draw.clear_flowers();
}

void tile_smap_data_t::render_scene_shadow_pass(point const &lpos) {
	terrain_tile_draw.draw_shadow_pass(lpos, tile);
}

bool tile_smap_data_t::needs_update(point const &lpos) {

	//return smap_data_t::needs_update(lpos);
	// FIXME: it would be better if we could just translate the shadow map when the scene shifts, but this seems fairly complex to track and get right
	int const new_dxoff(xoff - xoff2), new_dyoff(yoff - yoff2);
	bool const new_off(new_dxoff != dxoff || new_dyoff != dyoff);
	dxoff = new_dxoff; dyoff = new_dyoff;
	return (smap_data_t::needs_update(lpos) || new_off);
}


tile_t *get_tile_from_xy  (tile_xy_pair const &tp) {return terrain_tile_draw.get_tile_from_xy(tp);}
float update_tiled_terrain(float &min_camera_dist) {return terrain_tile_draw.update(min_camera_dist);}
void pre_draw_tiled_terrain() {terrain_tile_draw.pre_draw();}

colorRGBA get_inf_terrain_mod_color() {

	colorRGBA const colors[NUM_FIRE_MODES] = {WHITE, GREEN, RED, BLUE};
	assert(inf_terrain_fire_mode < NUM_FIRE_MODES);
	return colors[inf_terrain_fire_mode];
}

bool line_intersect_tiled_mesh_get_tile(point const &v1, point const &v2, point &p_int, tile_t *&tile) {

	float t(0.0);
	int xpos(0), ypos(0); // unused
	if (!terrain_tile_draw.line_intersect_mesh(v1, v2, t, tile, xpos, ypos)) return 0;
	p_int = v1 + t*(v2 - v1);
	return 1;
}

int get_tiled_terrain_tid_under_point(point const &pos) {
	return terrain_tile_draw.get_tid_under_point(pos);
}

void draw_tiled_terrain(bool reflection_pass) {

	render_models(0, model3d_offset.get_xlate()); // not sure where this goes
	//RESET_TIME;
	terrain_tile_draw.draw(reflection_pass);
	//glFinish(); PRINT_TIME("Tiled Terrain Draw"); //exit(0);

	if (inf_terrain_fire_mode != 0 && !reflection_pass) { // use a bool instead?
		point const v1(get_camera_pos()), v2(v1 + cview_dir*FAR_CLIP);
		point hit_pos;
		tile_t *tile(nullptr);

		if (line_intersect_tiled_mesh_get_tile(v1, v2, hit_pos, tile)) {
			draw_single_colored_sphere(hit_pos, 0.1, N_SPHERE_DIV, RED);

			// modification marker area rendered with a cylinder + stencil test to mask to mesh surface
			assert(tile != nullptr);
			float const tzmin(min(hit_pos.z, tile->get_zmin())), tzmax(max(hit_pos.z, tile->get_zmax())), dz(tzmax - tzmin);
			float const z1(tzmin - dz - SMALL_NUMBER), z2(tzmax + dz + SMALL_NUMBER);
			float const radius((cur_brush_param.get_radius() + 0.5)*HALF_DXY); // do we need a nonuniform x/y scale for non-square meshes?
			point const p1(hit_pos.x, hit_pos.y, z1);
			enable_blend();
			shader_t s;
			s.begin_color_only_shader(colorRGBA(get_inf_terrain_mod_color(), 0.2)); // make mostly transparent (wireframe?)
			glClear(GL_STENCIL_BUFFER_BIT);
			glDepthMask(GL_FALSE);
			glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			glEnable(GL_STENCIL_TEST);
			glStencilFunc(GL_ALWAYS, 0, ~0);
			glStencilOpSeparate(GL_FRONT, GL_INCR_WRAP, GL_INCR_WRAP, GL_KEEP);
			glStencilOpSeparate(GL_BACK,  GL_DECR_WRAP, GL_DECR_WRAP, GL_KEEP);
			draw_cylinder_at(p1, z2-z1, radius, radius, N_CYL_SIDES, 1); // with ends
			glStencilFunc(GL_NOTEQUAL, -1, ~0); // one front face only
			glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			draw_cylinder_at(p1, z2-z1, radius, radius, N_CYL_SIDES, 1); // with ends
			glDepthMask(GL_TRUE);
			glDisable(GL_STENCIL_TEST);
			disable_blend();
			s.end_shader();
		}
	}
	if (!reflection_pass && camera_surf_collide) {
		//int const tid(get_tiled_terrain_tid_under_point(get_camera_pos())); // TESTING
	}
}

void draw_tiled_terrain_lightning(bool reflection_pass) {terrain_tile_draw.update_lightning(reflection_pass);}
void clear_tiled_terrain() {terrain_tile_draw.clear();}
void reset_tiled_terrain_state() {terrain_tile_draw.clear_vbos_tids();}
void clear_tiled_terrain_shaders() {terrain_tile_draw.free_compute_shader();}
void draw_tiled_terrain_water(shader_t &s, float zval) {terrain_tile_draw.draw_water(s, zval);}
bool check_player_tiled_terrain_collision() {return terrain_tile_draw.check_player_collision();}
float get_tiled_terrain_water_level() {return (is_water_enabled() ? water_plane_z : terrain_tile_draw.get_actual_zmin());}


// *** heightmap modification and queries ***


bool line_intersect_tiled_mesh(point const &v1, point const &v2, point &p_int) {

	tile_t *tile(nullptr); // unused
	return line_intersect_tiled_mesh_get_tile(v1, v2, p_int, tile);
}


bool hmap_mod_enabled() {return (inf_terrain_fire_mode && using_tiled_terrain_hmap_tex());}

void change_inf_terrain_fire_mode(int val) {

	if (!using_tiled_terrain_hmap_tex()) return; // ignore
	inf_terrain_fire_mode = (inf_terrain_fire_mode + NUM_FIRE_MODES + val) % NUM_FIRE_MODES;
	string const modes[NUM_FIRE_MODES] = {"Look Only", "Increase Mesh Height", "Decrease Mesh Height", "Flatten Mesh"};
	print_text_onscreen(modes[inf_terrain_fire_mode], WHITE, 1.0, TICKS_PER_SECOND, 1); // 1 second
}

tile_t *get_tile_for_xy(int x, int y) {

	int const tsz(get_tile_size());
	if (x < 0) {x -= tsz-1;} // handle truncation toward lower integer
	if (y < 0) {y -= tsz-1;}
	return get_tile_from_xy(tile_xy_pair(x/tsz, y/tsz));
}

void inf_terrain_fire_weapon() {

	if (!hmap_mod_enabled()) return; // ignore
	//RESET_TIME;
	static float last_tfticks(0.0);
	if ((tfticks - last_tfticks) <= cur_brush_param.delay) return; // limit firing rate
	last_tfticks = tfticks;
	point const v1(get_camera_pos()), v2(v1 + cview_dir*FAR_CLIP);
	float t(0.0); // unused
	tile_t *tile(NULL);
	int xpos(0), ypos(0);
	if (!terrain_tile_draw.line_intersect_mesh(v1, v2, t, tile, xpos, ypos)) return;
	// Note: update is slow when trees are enabled
	unsigned shape(cur_brush_param.shape);
	if (inf_terrain_fire_mode == 3) {shape = ((shape == BSHAPE_CONST_SQ) ? BSHAPE_FLAT_SQ : BSHAPE_FLAT_CIR);} // enable a flattening shape
	float const delta_mag(cur_brush_param.get_delta_mag()*((inf_terrain_fire_mode == 1) ? 1.0 : -1.0));
	tex_mod_map_manager_t::hmap_val_t const base_delta(terrain_hmap_manager.scale_delta(delta_mag));
	tex_mod_map_manager_t::hmap_brush_t const brush(xpos, ypos, base_delta, cur_brush_param.get_radius(), shape);
	terrain_hmap_manager.apply_brush(brush, tile, 1); // cache
	//PRINT_TIME("Hmap Brush");
}

void inf_terrain_undo_hmap_mod() {

	if (!hmap_mod_enabled()) return;
	tex_mod_map_manager_t::hmap_brush_t brush;
	if (!terrain_hmap_manager.pop_last_brush(brush)) return;
	if (brush.is_flatten_brush()) return; // can't undo this brush since it's lossy
	// FIXME: won't work if clamping to min/max height occurred when applying the brush the first time
	brush.delta = -brush.delta; // invert
	terrain_hmap_manager.apply_brush(brush, get_tile_for_xy(brush.x, brush.y), 0); // don't cache
}


