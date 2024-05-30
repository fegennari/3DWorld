// 3D World - Tiled Landscape Mesh Generation/Drawing Code
// by Frank Gennari
// 9/26/10

#include "tiled_mesh.h"
#include "textures.h"
#include "shaders.h"
#include "openal_wrap.h"
#include "heightmap.h"
#include "profiler.h"


bool const DEBUG_TILES        = 0;
bool const DEBUG_TILE_BOUNDS  = 0;
bool const ENABLE_INST_PINE   = 1; // faster generation, lower GPU memory, slower rendering
bool const ENABLE_ANIMALS     = 1;
bool const USE_PARAMS_HSCALE  = 0;
bool const FLATTEN_BUILDING_TILE = 1; // removes terrain from the inside of buildings, but is slightly slower/higher memory usage and requires space between building and tile edge to prevent seams
int  const DITHER_NOISE_TEX   = NOISE_GEN_TEX;//PS_NOISE_TEX
unsigned const NORM_TEXELS    = 512;
unsigned const TILE_SMAP_START_TU_ID = 13;
float const FOG_DIST_TILES    = 1.45;
float const DRAW_DIST_TILES   = 1.5;
float const CREATE_DIST_TILES = 1.6;
float const CLEAR_DIST_TILES  = 1.6;
float const DELETE_DIST_TILES = 1.8;
float const GRASS_LOD_SCALE   = 15.0; // smaller = more grass detail
float const GRASS_DIST_SLOPE  = 0.25;
float const GRASS_THRESH      = 1.6;
float const SMAP_NEW_THRESH   = 1.2;
float const SMAP_DEL_THRESH   = 1.3;
float const SMAP_FADE_THRESH  = 1.5;
float const OCCLUDER_DIST     = 0.2;
float const FLOWER_REL_DIST   = 0.9; // flower view distance relative to grass view distance

int   const LIGHTNING_LIGHT = 2;
float const LIGHTNING_FREQ  = 200.0; // in ticks (1/40 s)
float const LITNING_TIME2   = 40.0;
float const LITNING_DIST    = 1.2;

unsigned const NUM_AO_DIRS  = 8; // Note: required to be 8 for adj tile calculation
unsigned const NUM_AO_STEPS = 8;
unsigned const AO_RAY_LEN(NUM_AO_STEPS*(NUM_AO_STEPS+1)/2); // 36

enum {FM_NONE, FM_INC_MESH, FM_DEC_MESH, FM_FLATTEN, FM_REM_TREES, FM_ADD_TREES, FM_REM_GRASS, FM_ADD_GRASS, NUM_FIRE_MODES};

struct clear_area_t : public cube_t {
	bool clear_next_frame;
	clear_area_t(cube_t const &region, bool cnf) : cube_t(region), clear_next_frame(cnf) {}
};


bool tt_lightning_enabled(0), check_tt_mesh_occlusion(1), shadow_maps_disabled(0);
unsigned inf_terrain_fire_mode(0); // none, increase height, decrease height
string read_hmap_modmap_fn, write_hmap_modmap_fn("heightmap.mod");
hmap_brush_param_t cur_brush_param;
tile_offset_t model3d_offset;
vector<clear_area_t> tile_smaps_to_clear;

extern bool inf_terrain_scenery, enable_tiled_mesh_ao, underwater, fog_enabled, volume_lighting, combined_gu, enable_depth_clamp, tt_triplanar_tex, use_grass_tess;
extern bool use_instanced_pine_trees, enable_tt_model_reflect, water_is_lava, tt_fire_button_down, flashlight_on, camera_in_building, rotate_trees;
extern unsigned grass_density, max_unique_trees, shadow_map_sz, erosion_iters_tt, num_rnd_grass_blocks, tiled_terrain_gen_heightmap_sz;
extern unsigned num_birds_per_tile, num_fish_per_tile, num_bflies_per_tile, room_geom_mem;
extern int DISABLE_WATER, display_mode, tree_mode, leaf_color_changed, ground_effects_level, animate2, iticks, num_trees, window_width, window_height, player_in_basement;
extern int invert_mh_image, is_cloudy, camera_surf_collide, show_fog, mesh_gen_mode, mesh_gen_shape, cloud_model, precip_mode, auto_time_adv, draw_model;
extern int player_in_elevator, player_in_attic;
extern float zmax, zmin, water_plane_z, mesh_scale, mesh_scale_z, vegetation, relh_adj_tex, grass_length, grass_width, fticks, cloud_height_offset, clouds_per_tile;
extern float ocean_wave_height, sm_tree_density, tree_density_thresh, atmosphere, cloud_cover, temperature, flower_density, FAR_CLIP, biome_x_offset;
extern float smap_thresh_scale, tt_grass_scale_factor, pond_max_depth;
extern double tfticks;
extern point sun_pos, moon_pos, surface_pos;
extern vector3d wind;
extern cube_t grass_exclude1, grass_exclude2, building_occluder;
extern water_params_t water_params;
extern char *mh_filename_tt;
extern float h_dirt[];
extern tree_data_manager_t tree_data_manager;
extern pt_line_drawer tree_scenery_pld;
extern tree_placer_t tree_placer;

bool enable_terrain_env(ENABLE_TERRAIN_ENV);
void set_water_plane_uniforms(shader_t &s);
void create_pine_tree_instances();
unsigned get_tree_inst_gpu_mem();
void setup_detail_normal_map(shader_t &s, float tscale);
void draw_distant_mesh_bottom(float terrain_zmin);
bool no_grass_under_buildings();
bool check_buildings_no_grass(point const &pos);
colorRGBA get_avg_color_for_landscape_tex(unsigned id); // defined later in this file
void building_gameplay_action_key(int mode, bool mouse_wheel);
bool player_cant_see_outside_building();
bool check_cube_occluded(cube_t const &cube, vect_cube_t const &occluders, point const &viewer);
void get_city_grass_coll_cubes(cube_t const &region, vect_cube_t &out, vect_cube_t &out_bt);
int check_city_contains_overlaps(cube_t const &query);
bool check_inside_city(point const &pos, float radius);
cube_t get_city_bcube_overlapping(cube_t const &c);
void show_gpu_mem_info();


float get_inf_terrain_fog_dist() {return FOG_DIST_TILES*get_scaled_tile_radius();}
float get_draw_tile_dist  () {return DRAW_DIST_TILES*get_scaled_tile_radius();}
float get_grass_thresh    () {return GRASS_THRESH*tt_grass_scale_factor*get_tile_width();}
float get_grass_blend_dist() {return tt_grass_scale_factor/GRASS_DIST_SLOPE;}
float get_grass_thresh_pad() {return (get_grass_thresh() + get_grass_blend_dist());}
bool is_water_enabled     () {return (!DISABLE_WATER && (display_mode & 0x04) != 0);}
bool pine_trees_enabled   () {return (((tree_mode & 2) && vegetation > 0.0) || !tree_placer.sm_blocks.empty());} // and palm trees
bool decid_trees_enabled  () {return (((tree_mode & 1) && vegetation > 0.0) || !tree_placer.   blocks.empty());}
bool any_trees_enabled    () {return (pine_trees_enabled() || decid_trees_enabled());}
bool scenery_enabled      () {return (inf_terrain_scenery && SCENERY_THRESH > 0.0);}
bool gen_grass_map        () {return (GRASS_THRESH > 0.0 && grass_density > 0 && vegetation > 0.0);}
bool is_grass_enabled     () {return ((display_mode & 0x02) && gen_grass_map());}
bool clouds_enabled       () {return ((display_mode & 0x40) == 0 && atmosphere > 0.0);}
bool cloud_shadows_enabled() {return (ground_effects_level >= 2 && clouds_enabled());}
bool mesh_shadows_enabled () {return (ground_effects_level >= 1);}
bool is_distance_mode     () {return ((display_mode & 0x10) != 0);}
bool nonunif_fog_enabled  () {return (show_fog && is_distance_mode());}
bool enable_ocean_waves   () {return ((display_mode & 0x0100) != 0 && wind.mag() > TOLERANCE);}
bool draw_distant_water   () {return (is_water_enabled() && is_distance_mode() && far_clip_ratio > 1.1);}
bool enable_tree_dlights  () {return ((is_night() && have_cities()) || flashlight_on);} // enable for city night lights
float get_tt_fog_top      () {return (nonunif_fog_enabled() ? (zmax + (zmax - zmin)) : (zmax + FAR_CLIP));}
float get_tt_fog_bot      () {return (nonunif_fog_enabled() ? zmax : (zmax + FAR_CLIP));}
float get_tt_cloud_level  () {return 0.5f*(get_tt_fog_bot() + get_tt_fog_top());}
float get_smap_atten_val  () {return SMAP_FADE_THRESH*smap_thresh_scale*get_tile_width();}
float get_tile_smap_dist  () {return get_smap_atten_val();}
float get_max_sea_level   () {return (get_water_z_height() + ocean_wave_height);}
unsigned get_tile_size    () {return MESH_X_SIZE;}

bool use_water_plane_tess () {
	if (!enable_ocean_waves() || cloud_model != 0 || draw_distant_water()) return 0; // hack to use cloud_model (F10)
	static bool tess_enabled(1);
	if (tess_enabled && !check_for_tess_shader()) {tess_enabled = 0;} // disable tess - not supported
	return tess_enabled;
}

vector3d get_tiled_terrain_model_xlate() {return model3d_offset.get_xlate();}
vector3d get_camera_coord_space_xlate () {return vector3d((world_mode == WMODE_INF_TERRAIN) ? get_tiled_terrain_model_xlate() : zero_vector);}

bool enable_instanced_pine_trees() {
	if (use_instanced_pine_trees) return 1;
	if (have_cities())            return 0; // disable if there are cities because then we can't create custom pine trees for cities
	float const ntrees_mult(vegetation*sm_tree_density*tree_density_thresh*tree_scale*tree_scale);
	return (ENABLE_INST_PINE && (tree_mode & 2) && ntrees_mult >= ((tree_mode == 3) ? 3 : 4) && max_unique_trees > 0); // enable when there are lots of pine/palm trees
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
		dfavg = 0.5f*(zf + z0) - zm;
	}
	float const fog_dist((zc - zm)*sqrt(max((FD*FD/(dfavg*dfavg) - 1.0), 0.0)));
	//cout << "FD: " << FD << ", fog_dist: " << fog_dist << ", FAR_CLIP: " << FAR_CLIP << ", final: " << max(FAR_CLIP, 1.1f*fog_dist) << endl;
	return max(uniform_fog_far_clip, 1.1f*fog_dist); // add 10% padding
}


grass_tile_manager_t grass_tile_manager;

void update_tiled_terrain_grass_vbos() {grass_tile_manager.clear_vbo();}
void update_tiled_grass_colors() {grass_tile_manager.clear();} // regenerate grass


#define BILINEAR_INTERP(arr, var, x, y) (y*(x*arr[1][1].var + (1.0f-x)*arr[1][0].var) + (1.0f-y)*(x*arr[0][1].var + (1.0f-x)*arr[0][0].var))


// *** heightmap management ***


class tiled_terrain_hmap_manager_t : public terrain_hmap_manager_t {

	tile_t *cur_tile=nullptr;
	bool modified[3][3];

public:
	tiled_terrain_hmap_manager_t() {clear_modified();}
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
	void flatten_region(cube_t const &cube) {
		// Note: to be applied before tiles are generated so that they don't need to be invalidated
		// Note: assumes unscaled mesh (mesh_scale == 1)
		int const x1(floor((cube.x1() + X_SCENE_SIZE)*DX_VAL_INV)), y1(floor((cube.y1() + Y_SCENE_SIZE)*DY_VAL_INV));
		int const x2(ceil ((cube.x2() + X_SCENE_SIZE)*DX_VAL_INV)), y2(ceil ((cube.y2() + Y_SCENE_SIZE)*DY_VAL_INV));
		int cx1(x1), cy1(y1), cx2(x2+1), cy2(y2+1); // Note: cx2/cy2 are one past the end; this is needed for proper mirror clamping and empty range early termination optimization
		bool const allow_wrap(0); // while this mostly works with buildings, it ruins roads and can cause mesh seams
		if (!clamp_xy(cx1, cy1, 0.0, 0.0, allow_wrap) || !clamp_xy(cx2, cy2, 0.0, 0.0, allow_wrap)) return; // off the texture, skip
		
		if (allow_wrap) { // handle swapping due to X/Y mirroring
			if (cx2 < cx1) {swap(cx1, cx2);}
			if (cy2 < cy1) {swap(cy1, cy2);}
		}
		else {assert(cx1 >= 0 && cy1 >= 0 && cx1 <= cx2 && cy1 <= cy2);}
		if (cx1 == cx2 || cy1 == cy2) return; // empty range optimization
		point const center(cube.get_cube_center());
		float xc((center.x + X_SCENE_SIZE)*DX_VAL_INV + 0.5), yc((center.y + Y_SCENE_SIZE)*DY_VAL_INV + 0.5); // convert from real to index space
		int const xlo(floor(xc)), ylo(floor(yc)), xhi(ceil(xc)), yhi(ceil(yc));
		float const xv(xc - xlo), yv(yc - ylo); // use cubic_interpolate()?
		int const height_val(yv*(xv*get_clamped_pixel_value(xhi, yhi, 0) + (1.0f-xv)*get_clamped_pixel_value(xlo, yhi, 0)) +
			          (1.0f-yv)*(xv*get_clamped_pixel_value(xhi, ylo, 0) + (1.0f-xv)*get_clamped_pixel_value(xlo, ylo, 0))); // linear interpolation
		tex_mod_map_manager_t::mod_elem_t elem(cx1, cy1, height_val);

		for (elem.y = cy1; elem.y < cy2; ++elem.y) {
			for (elem.x = cx1; elem.x < cx2; ++elem.x) {modify_height(elem, 0);} // not wrapped
		}
	}
	virtual bool modify_height_value(int x, int y, hmap_val_t val, bool is_delta, float fract_x, float fract_y, bool allow_wrap=1) {
		int clamped_x(x), clamped_y(y);
		if (!clamp_xy(clamped_x, clamped_y, fract_x, fract_y, allow_wrap)) return 0;
		assert(clamped_x >= 0 && clamped_y >= 0);
		modify_height(tex_mod_map_manager_t::mod_elem_t(clamped_x, clamped_y, val), is_delta); // Note: *not* cached at this level
		if (cur_tile) {cur_tile->fill_adj_mask(modified, x, y);}
		return 1;
	}
};


tiled_terrain_hmap_manager_t terrain_hmap_manager;


bool using_tiled_terrain_hmap_tex() {return (world_mode == WMODE_INF_TERRAIN && terrain_hmap_manager.enabled());}
bool using_hmap_with_detail      () {return (using_tiled_terrain_hmap_tex() && mesh_scale < 0.75);}

float get_tiled_terrain_height_tex(float xval, float yval, bool nearest_texel) {
	return (nearest_texel ? terrain_hmap_manager.get_nearest_height(xval, yval) : terrain_hmap_manager.interpolate_height(xval, yval));
}
vector3d get_tiled_terrain_height_tex_norm(int x, int y) {return terrain_hmap_manager.get_norm(x, y);}

bool read_default_hmap_modmap() {

	if (read_hmap_modmap_fn.empty()) return 0;
	if (!terrain_hmap_manager.read_and_apply_mod(read_hmap_modmap_fn)) return 0;
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

tile_t::tile_t() : decid_trees(tree_data_manager) {}

tile_t::tile_t(unsigned size_, int x, int y) : size(size_), stride(size+1), zvsize(stride+1),
	deltax(DX_VAL), deltay(DY_VAL), mesh_off(xoff-xoff2, yoff-yoff2), decid_trees(tree_data_manager)
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

float tile_t::get_draw_priority() const {
	return (p2p_dist_xy(get_camera_pos(), get_center()) + (is_visible() ? 0.0 : FAR_CLIP)); // prioritize visible tiles
}


void tile_t::update_terrain_params() { // setup biomes

	float const dirt_mult(1.0), veg_mult(5.0);
	float const xv1(get_xval(x1)), xv2(xv1 + (x2-x1)*deltax), yv1(get_yval(y1)), yv2(yv1 + (y2-y1)*deltay);

	for (unsigned yp = 0; yp < 2; ++yp) { // sample at 4 tile corners and interpolate across the tile
		for (unsigned xp = 0; xp < 2; ++xp) {
			terrain_params_t &param(params[yp][xp]);
			float const xv(mesh_scale*(xp ? xv2 : xv1) + biome_x_offset), yv(mesh_scale*(yp ? yv2 : yv1));

			if (USE_PARAMS_HSCALE) {
				param.hoff   = eval_mesh_sin_terms(0.4*xv+123, 0.4*yv+456);
				param.hscale = min(2.0f, max(0.5f, 0.5f*fabs(eval_mesh_sin_terms(0.8*xv+789, 0.8*yv+111))));
			}
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

	unsigned mem(pine_trees.get_gpu_mem() + decid_trees.get_gpu_mem() + scenery.get_gpu_mem() + flowers.get_gpu_mem());
	unsigned const num_texels(stride*stride);
	if (weight_tid > 0) {mem += 4*num_texels;} // 4 bytes per texel (RGBA8)
	if (height_tid > 0) {mem += 4*num_texels;} // 4 bytes per texel (F32)
	if (normal_tid > 0) {mem += 4*num_texels;} // 4 bytes per texel (RGBA8)
	if (shadow_tid > 0) {mem += 4*num_texels;} // 4 bytes per texel (RGBA8)
	mem += get_smap_mem(); // FBO textures
	return mem;
}

unsigned tile_t::get_smap_mem() const {
	unsigned mem(0);
	for (unsigned i = 0; i < smap_data.size(); ++i) {mem += smap_data[i].get_gpu_mem();}
	return mem;
}

unsigned tile_t::count_shadow_maps() const {
	unsigned num(0);
	for (unsigned i = 0; i < smap_data.size(); ++i) {num += smap_data[i].is_allocated();}
	return num;
}


void tile_t::clear() {

	clear_vbo_tid(nullptr);
	tree_map.clear();
	mesh_weight_data.clear();
	weight_data.clear();
	zvals.clear();
	clear_shadows();
	pine_trees.clear_all();
	decid_trees.clear();
	scenery.clear();
	grass_blocks.clear();
	clear_flowers();
}

void tile_t::clear_shadows(bool clear_sun, bool clear_moon, bool no_clear_adj) {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if ((l == LIGHT_SUN && !clear_sun) || (l == LIGHT_MOON && !clear_moon)) continue;
		smask[l].clear();
		
		if (!no_clear_adj) { // this is not necessary, but is an optimization for low sun angles to avoid recomputation of shadows on a tile when its neighbors have changed
			for (unsigned d = 0; d < 2; ++d) {sh_out[l][d].clear();}
		}
	}
	sun_shadows_invalid  = clear_sun;
	moon_shadows_invalid = clear_moon;
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
	free_texture(normal_tid);
	free_texture(shadow_tid);
	gen_tsize = 0;
}


bool setup_height_gen(mesh_xy_grid_cache_t &height_gen, float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values, bool no_wait=0) {

	bool const add_detail(using_hmap_with_detail());
	if (!add_detail && using_tiled_terrain_hmap_tex()) return 1; // nothing to do
	float const xy_scale(add_detail ? HMAP_DETAIL_SCALE : 1.0);
	bool const results_avail(height_gen.build_arrays(xy_scale*x0, xy_scale*y0, xy_scale*dx, xy_scale*dy, nx, ny, cache_values, 0, no_wait));
	height_gen.enable_glaciate();
	return results_avail;
}


bool tile_t::create_zvals(mesh_xy_grid_cache_t &height_gen, bool no_wait) {

	//timer_t timer("Create Zvals");
	inside_city = check_city_contains_overlaps(get_mesh_bcube_global());
	if (enable_terrain_env) {update_terrain_params();}
	zvals.resize(zvsize*zvsize);
	mzmin =  FAR_DISTANCE;
	mzmax = -FAR_DISTANCE;
	unsigned const block_size(zvsize/4), context_sz(stride + 2*AO_RAY_LEN);
	bool const using_hmap(using_tiled_terrain_hmap_tex()), add_detail(using_hmap_with_detail()); // add procedural detail to heightmap

	// When using AO + GPU noise generation, it's faster to compute the AO + context and clip the zvals from this rather than making two separate compute calls (one without blocking)
	if (enable_tiled_mesh_ao && !using_hmap && mesh_gen_mode >= MGEN_SIMPLEX_GPU) {
		bool results_ready(setup_height_gen(height_gen, get_xval(x1 - AO_RAY_LEN), get_yval(y1 - AO_RAY_LEN), deltax, deltay, context_sz, context_sz, 0, no_wait)); // cache_values=0
		if (!results_ready) {assert(no_wait); return 0;} // cached heights are not yet ready
		ao_zvals.resize(context_sz*context_sz);

#pragma omp parallel for schedule(static,1)
		for (int y = 0; y < (int)context_sz; ++y) {
			for (unsigned x = 0; x < context_sz; ++x) {ao_zvals[y*context_sz + x] = height_gen.eval_index(x, y);}
		}
	}
	else {
		bool results_ready(setup_height_gen(height_gen, get_xval(x1), get_yval(y1), deltax, deltay, zvsize, zvsize, 0, no_wait)); // cache_values=0
		if (!results_ready) {assert(no_wait); return 0;} // cached heights are not yet ready
	}
	float const xy_mult(1.0/float(size)), wpz_max(get_max_sea_level());

#pragma omp parallel for schedule(static,1)
	for (int y = 0; y < (int)zvsize; ++y) {
		for (unsigned x = 0; x < zvsize; ++x) {
			float &zval(zvals[y*zvsize + x]);

			if (using_hmap) {
				zval = terrain_hmap_manager.get_clamped_height((x1 + x), (y1 + y));
				if (add_detail) {zval += HMAP_DETAIL_MAG*height_gen.eval_index(x, y);} // less hard-coded - scale by delta between adjacent zvals?
			}
			else {
				if (!ao_zvals.empty()) {zval = ao_zvals[(y + AO_RAY_LEN)*context_sz + (x + AO_RAY_LEN)];} // use AO zvals
				else                   {zval = height_gen.eval_index(x, y);} // use height gen

				if (USE_PARAMS_HSCALE) {
					float const xv(float(x)*xy_mult), yv(float(y)*xy_mult);
					zval = BILINEAR_INTERP(params, hoff, xv, yv) + BILINEAR_INTERP(params, hscale, xv, yv)*zval;
				}
			}
		} // for x
	} // for y
	if (!using_hmap) {apply_erosion(&zvals.front(), zvsize, zvsize, zmin, erosion_iters_tt);} // heightmap is eroded during load

	for (unsigned yy = 0; yy < 4; ++yy) {
		for (unsigned xx = 0; xx < 4; ++xx) {
			unsigned const x_end((xx+1)*block_size), y_end((yy+1)*block_size); // last row/column is skipped because it's not rendered
			assert(x_end < zvsize && y_end < zvsize);
			float &szmin(sub_zmin[yy][xx]), &szmax(sub_zmax[yy][xx]);
			szmin = FAR_DISTANCE; szmax = -FAR_DISTANCE;

			for (unsigned y = yy*block_size; y <= y_end; ++y) {
				for (unsigned x = xx*block_size; x <= x_end; ++x) {
					float const z(zvals[y*zvsize + x]);
					szmin = min(szmin, z); szmax = max(szmax, z);

					if (z < wpz_max) { // can be underwater - update water bbox
						wx1 = min(wx1, x1+int(x)); wy1 = min(wy1, y1+int(y));
						wx2 = max(wx2, x1+int(x)); wy2 = max(wy2, y1+int(y));
					}
				} // for x
			} // for y
			max_eq(mesh_dz, (szmax - szmin));
			mzmin = min(mzmin, szmin);
			mzmax = max(mzmax, szmax);
		} // for xx
	} // for yy
	assert(mzmin <= mzmax);
	radius = 0.5*sqrt((deltax*deltax + deltay*deltay)*size*size + (mzmax - mzmin)*(mzmax - mzmin));
	ptzmax = dtzmax = mzmin; // no trees yet
	if (!can_have_trees()) {no_trees = 1;} // mark as no_trees so that trees don't pop when water is disabled later
	if (DEBUG_TILES) {cout << "new tile coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;}
	return 1; // results are ready
}

void tile_t::get_z_minmax_for_area(point const &pos, float radius, float &zmin, float &zmax) const {

	float const rx1(pos.x - radius), ry1(pos.y - radius), rx2(pos.x + radius), ry2(pos.y + radius);
	unsigned const ix1(max(0, (get_xpos_round_down(rx1) - x1))); // stride = zvsize-1
	unsigned const iy1(max(0, (get_ypos_round_down(ry1) - y1)));
	unsigned const ix2(min(stride, (unsigned)(get_xpos_round_down(rx2) - x1 + 1)));
	unsigned const iy2(min(stride, (unsigned)(get_ypos_round_down(ry2) - y1 + 1)));
	assert(ix1 <= ix2 && iy1 <= iy2); // must be inside the tile

	for (unsigned y = iy1; y <= iy2; ++y) {
		for (unsigned x = ix1; x <= ix2; ++x) {
			float const z(zvals[y*zvsize + x]);
			min_eq(zmin, z);
			max_eq(zmax, z);
		}
	}
}

void get_tile_z_minmax_for_area(tile_t const &tile, point const &pos, float radius, float &zmin, float &zmax) {
	tile.get_z_minmax_for_area(pos, radius, zmin, zmax);
}

float tile_t::get_zval_at(float x, float y, bool in_global_space) const {

	assert(!zvals.empty());
	cube_t bcube(in_global_space ? get_mesh_bcube_global() : get_mesh_bcube());
	float const xx((x - bcube.d[0][0])/(bcube.d[0][1] - bcube.d[0][0])), yy((y - bcube.d[1][0])/(bcube.d[1][1] - bcube.d[1][0]));
	float const xp((x2 - x1)*CLIP_TO_01(xx)), yp((y2 - y1)*CLIP_TO_01(yy));
	int const x0((int)xp), y0((int)yp);
	float const xpi(xp - (float)x0), ypi(yp - (float)y0); // always positive
	assert(x0 >= 0 && y0 >= 0 && x0+1 < (int)zvsize && y0+1 < (int)zvsize); //return zvals[y*zvsize + x];
	unsigned const ix(y0*zvsize + x0);
	return (1.0f - xpi)*((1.0f - ypi)*zvals[ix] + ypi*zvals[ix + zvsize]) + xpi*((1.0f - ypi)*zvals[ix + 1] + ypi*zvals[ix + zvsize + 1]);
}


// *** shadows + AO lighting ***

void tile_t::calc_mesh_ao_lighting() {

	//timer_t timer("Calc Tile AO Lighting");
	// caclulate ray step directions
	tile_xy_pair ao_dirs[NUM_AO_DIRS]; // 0  1  2  3  4  5  6  7
	unsigned ix(0);

	for (int y = -1; y <= 1; ++y) {
		for (int x = -1; x <= 1; ++x) {
			if (x != 0 || y != 0) {ao_dirs[ix++] = tile_xy_pair(x, y);}
		}
	}
	assert(ix == NUM_AO_DIRS);
	assert(AO_RAY_LEN <= size);

	// create context zvals, which may overlap with other tiles (that need not be created at this point)
	unsigned const context_sz(stride + 2*AO_RAY_LEN);
	bool const using_hmap(using_tiled_terrain_hmap_tex()), add_detail(using_hmap_with_detail()), use_ao_zvals(!ao_zvals.empty());
	vector<float> czv;
	mesh_xy_grid_cache_t height_gen;
	
	if (use_ao_zvals) {czv.swap(ao_zvals);} // use precomputed values, will clear ao_zvals at the end
	else {
		czv.resize(context_sz*context_sz);
		setup_height_gen(height_gen, get_xval(x1 - AO_RAY_LEN), get_yval(y1 - AO_RAY_LEN), deltax, deltay, context_sz, context_sz, 0); // cache_values=0
	}
	float const dz(0.5*HALF_DXY);
	ao_lighting.resize(stride*stride);

#pragma omp parallel
	{
		if (!use_ao_zvals) {
#pragma omp for schedule(static,1)
			for (int y = 0; y < (int)context_sz; ++y) {
				for (int x = 0; x < (int)context_sz; ++x) {
					int const xv(x - AO_RAY_LEN), yv(y - AO_RAY_LEN);
					float &zv(czv[y*context_sz + x]);
					if (xv >= 0 && yv >= 0 && xv < (int)zvsize && yv < (int)zvsize) {zv = zvals[yv*zvsize + xv];}
					else if (using_hmap) {
						zv = terrain_hmap_manager.get_clamped_height((x1 + xv), (y1 + yv));
						if (add_detail) {zv += HMAP_DETAIL_MAG*height_gen.eval_index(x, y);}
					}
					else {zv = height_gen.eval_index(x, y);} // Note: not using hoff/hscale here since they are undefined outside the tile bounds
				}
			}
		}
		// calculate ao_lighting values by casting rays through the mesh zvals
#pragma omp for schedule(static,1)
		for (int y = 0; y < (int)stride; ++y) {
			for (int x = 0; x < (int)stride; ++x) {
				unsigned atten(0);

				for (unsigned d = 0; d < NUM_AO_DIRS; ++d) {
					float z0(zvals[y*zvsize + x]);
					tile_xy_pair step(ao_dirs[d]);
					tile_xy_pair v(x, y);

					for (unsigned s = 0; s < NUM_AO_STEPS; ++s) {
						v    += step;
						z0   += dz;
						//step += step; // multiply by 2 for exponential step size
						step += ao_dirs[d]; // linear increase (Note: must agree with max_ray_length)
						int const xv(v.x + AO_RAY_LEN), yv(v.y + AO_RAY_LEN);
						//assert(xv >= 0 && yv >= 0 && xv < (int)context_sz && yv < (int)context_sz);
						
						if (czv[yv*context_sz + xv] > z0) { // hit a higher point
							atten += (NUM_AO_STEPS - s); // Note: ambient obscurance - uses actual distance to occluder
							break;
						}
					} // for s
				} // for d
				assert(atten <= NUM_AO_DIRS*NUM_AO_STEPS);
				float const ao_scale(1.0 - float(atten)/float(NUM_AO_DIRS*NUM_AO_STEPS));
				ao_lighting[y*stride + x] = (unsigned char)(255.0*ao_scale);
			} // for x
		} // for y
	}
}


void tile_t::calc_shadows_for_light(unsigned l) {

	if (is_distant) return; // Note: can be made to work, but won't work as-is
	assert(!smask[l].empty());
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
	((l == LIGHT_SUN) ? sun_shadows_invalid : moon_shadows_invalid) = 1;
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
	} // end while()
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
	int const rval(max(int(tradius/deltax), int(tradius/deltay)) + 1), rval_sq(rval*rval);
	int const x1(max(0, xc-rval)), y1(max(0, yc-rval)), x2(min((int)size, xc+rval)), y2(min((int)size, yc+rval));
	float const scale(0.6/rval);
	bool updated(0);

	for (int y = y1; y <= y2; ++y) {
		for (int x = x1; x <= x2; ++x) {
			float const dx(abs(x - xc)), dy(abs(y - yc)), dist_sq(dx*dx + dy*dy);
			if (dist_sq > rval_sq) continue;
			float const mult(0.2 + 0.8*scale*sqrt(dist_sq));
			tree_map_val &val(tree_map[y*stride + x]);
			val.ao *= mult;
			val.sh *= mult;
			updated = 1;
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
		sun_shadows_invalid = moon_shadows_invalid = 1;
		if (has_any_grass) {recalc_tree_grass_weights = 1;} // Note: may be slow, and doesn't have a big impact
	}
}

template<typename T> void tile_t::apply_ao_shadows_for_tree_group(T const &trees, tile_offset_t const &toff, bool no_adj_test, float rscale) {

	point const pt_off(toff.subtract_from(mesh_off));
	cube_t const bcube(get_mesh_bcube());

	for (typename T::const_iterator i = trees.begin(); i != trees.end(); ++i) {
		point const pt(i->get_center() + pt_off);
		float const tr(rscale*i->get_radius());
		if (no_adj_test && (pt.x+tr < bcube.d[0][0] || pt.x-tr > bcube.d[0][1] || pt.y+tr < bcube.d[1][0] || pt.y-tr > bcube.d[1][1])) continue;
		add_tree_ao_shadow(pt, tr, no_adj_test);
	}
}

void tile_t::apply_ao_shadows_for_trees(tile_t const *const tile, bool no_adj_test) {

	if (can_have_pine_palm_trees()) {
		assert(pine_trees_generated());
		apply_ao_shadows_for_tree_group(tile->pine_trees, tile->ptree_off, no_adj_test, 1.8);
	}
	if (can_have_decid_trees()) {
		assert(decid_trees.was_generated());
		apply_ao_shadows_for_tree_group(tile->decid_trees, tile->dtree_off, no_adj_test, 0.5);
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
	//timer_t timer("Tree AO Shadows");
	tree_map.resize(0);
	tree_map.resize(stride*stride);
	bool const no_adj_test(trmax < min(deltax, deltay));
	apply_ao_shadows_for_trees(this, no_adj_test);
}


void tile_t::check_shadow_map_and_normal_texture(bool no_push) {

	if (!normal_tid) {
		setup_texture(normal_tid, 0, 0, 0, 0, 0);
		upload_normal_texture(0); // created once, never updated (so never valid here)
	}
	bool const tid_is_valid(shadow_tid != 0);
	bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6), mesh_shadows(mesh_shadows_enabled());
	bool const update_sun (has_sun  && (sun_shadows_invalid  || smask[LIGHT_SUN ].empty()));
	bool const update_moon(has_moon && (moon_shadows_invalid || smask[LIGHT_MOON].empty()));
	if (tid_is_valid && !(update_sun || update_moon)) return; // up-to-date
	//timer_t timer("Shadow Map Texture Update");
	if (!tid_is_valid) {setup_texture(shadow_tid, 0, 0, 0, 0, 0); sun_shadows_invalid = moon_shadows_invalid = 1;}
	assert(has_sun || has_moon);
	if (mesh_shadows) {calc_shadows(update_sun, update_moon, no_push);}
	if (enable_tiled_mesh_ao && ao_lighting.empty()) {calc_mesh_ao_lighting();}
	upload_shadow_map_texture(tid_is_valid);
	sun_shadows_invalid = moon_shadows_invalid = 0;
}


// Note: all of these textures are really RGB, but we upload them as RGBA for proper 4-byte alignment (since they are a power of 2 + 1)
void create_or_update_texture(unsigned &tid, bool tid_is_valid, unsigned stride, vector<unsigned char> const &data) {

	bind_2d_texture(tid);

	if (tid_is_valid) { // overwrite old data
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, stride, stride, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
	}
	else { // allocate and write
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, stride, stride, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data.front());
	}
}

void tile_t::upload_normal_texture(bool tid_is_valid) {

	//timer_t timer("Create Normal Texture");
	vector<unsigned char> normal_data(4*stride*stride, 0);
	min_normal_z = 1.0;

	for (unsigned y = 0; y < stride; ++y) {
		for (unsigned x = 0; x < stride; ++x) {
			unsigned const ix(y*stride + x), ix2(y*zvsize + x), ix_off(4*ix);
			vector3d const norm(get_norm(ix2));
			min_normal_z = min(min_normal_z, norm.z);
			UNROLL_3X(normal_data[ix_off+i_] = (unsigned char)(127.0*(norm[i_] + 1.0)););
		}
	}
	create_or_update_texture(normal_tid, tid_is_valid, stride, normal_data);
}

void tile_t::upload_shadow_map_texture(bool tid_is_valid) {

	//timer_t timer("Create Shadow Map Texture");
	bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6), mesh_shadows(mesh_shadows_enabled());
	vector<unsigned char> shadow_data(4*stride*stride, 0);
	vector<unsigned char> const &cur_smask(smask[has_sun ? (unsigned)LIGHT_SUN : (unsigned)LIGHT_MOON]);
	if (mesh_shadows && has_sun ) {assert(!smask[LIGHT_SUN ].empty());}
	if (mesh_shadows && has_moon) {assert(!smask[LIGHT_MOON].empty());}
	float const lfs(5.0*(light_factor - 0.4)); // 0: all moon, 1: all sun

	for (unsigned y = 0; y < stride; ++y) { // Note: shadow texture is stored as {mesh_shadow, tree_shadow, ambient_occlusion}
		for (unsigned x = 0; x < stride; ++x) {
			unsigned const ix(y*stride + x), ix2(y*zvsize + x), ix_off(4*ix);
			// 67% ambient if AO lighting is disabled (to cancel out with the scale by 1.5 in the shaders)
			unsigned char const base_ao(ao_lighting.empty() ? 170 : ao_lighting[ix]); // Note: must always do this so that adj tiles can update tree AO at tile borders
			if (tree_map.empty() || tree_map[ix].ao == 255) {shadow_data[ix_off+2] = base_ao;} // no tree AO
			else {shadow_data[ix_off+2] = (unsigned char)(base_ao * (0.3f + 0.7f*tree_map[ix].ao/255.0f));}
			unsigned char shadow_val(255);

			if (!mesh_shadows) {} // do nothing
			else if (has_sun && has_moon) {
				bool const sun_en( (smask[LIGHT_SUN ][ix2] & SHADOWED_ALL) == 0);
				bool const moon_en((smask[LIGHT_MOON][ix2] & SHADOWED_ALL) == 0);
				shadow_val *= (lfs*sun_en + (1.0f - lfs)*moon_en);
			}
			else if (cur_smask[ix2] & SHADOWED_ALL) {shadow_val = 0;} // fully in shadow
			shadow_data[ix_off+0] = shadow_val; // mesh shadow
			shadow_data[ix_off+1] = (tree_map.empty() ? 255 : (unsigned char)(63.75f + 0.75f*tree_map[ix].sh)); // 75% tree shadow: fully lit if not nearby trees
		}
	}
	create_or_update_texture(shadow_tid, tid_is_valid, stride, shadow_data);
}

tile_smap_data_t tile_shadow_map_manager::new_smap_data(unsigned tu_id, tile_t *tile, unsigned light, unsigned lod_level) {
	assert(tile != nullptr);
	assert(light < NUM_LIGHT_SRC);
	assert(lod_level < NUM_SMAP_LODS);
	unsigned const target_smap_sz(shadow_map_sz >> lod_level);
	auto &cur_free_list(free_list[light][lod_level]);
	if (cur_free_list.empty()) {return tile_smap_data_t(tu_id, target_smap_sz, lod_level, tile);}
	smap_data_state_t const state(cur_free_list.back());
	cur_free_list.pop_back();
	return tile_smap_data_t(tu_id, target_smap_sz, lod_level, tile, state);
}

void tile_shadow_map_manager::release_smap_data(tile_smap_data_t &smd, unsigned light) {
	if (!smd.is_allocated()) return;
	assert(light < NUM_LIGHT_SRC);
	assert(smd.lod_level < NUM_SMAP_LODS);
	free_list[light][smd.lod_level].push_back(smd); // adds the base class, which contains the GL state
	smd.disown(); // no longer owned by smd
}

void tile_shadow_map_manager::clear_context() {
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		for (unsigned L = 0; L < NUM_SMAP_LODS; ++L) {
			auto &cur_free_list(free_list[l][L]);
			for (auto i = cur_free_list.begin(); i != cur_free_list.end(); ++i) {i->free_gl_state();}
			cur_free_list.clear();
		}
	}
}

unsigned get_smap_bytes_per_pixel();

unsigned tile_shadow_map_manager::get_free_list_mem_usage() const {
	unsigned mem(0);

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		for (unsigned L = 0; L < NUM_SMAP_LODS; ++L) {
			unsigned const tex_size(shadow_map_sz >> L);
			mem += get_smap_bytes_per_pixel()*tex_size*tex_size*free_list[l][L].size();
		}
	}
	return mem;
}

cube_t tile_t::get_shadow_bcube() const {
	vector3d const b_ext(get_buildings_max_extent()); // what about bridges overlapping this tile?
	vector2d const road_len(get_road_max_len());
	float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2));
	float const x_ext(max(max(0.5f*road_len.x, b_ext.x), trmax)), y_ext(max(max(0.5f*road_len.y, b_ext.y), trmax));
	return cube_t(xv1-x_ext, xv1+(x2-x1)*deltax+x_ext, yv1-y_ext, yv1+(y2-y1)*deltay+y_ext, mzmin-BCUBE_ZTOLER, max(get_tile_zmax()+BCUBE_ZTOLER, mzmax+b_ext.z));
}

void tile_t::draw_smap_debug_vis(shader_t &s) const {
	if (smap_data.empty()) return; // no active shadow maps
	colorRGBA const lod_colors[6] = {RED, ORANGE, YELLOW, GREEN, BLUE, PURPLE}; // rainbow
	s.set_cur_color(lod_colors[min(smap_lod_level, 5U)]);
	draw_simple_cube(get_shadow_bcube(), 0);
	cout << (shadow_map_sz >> smap_lod_level) << " ";
}

unsigned calc_max_smap_lod() {
	unsigned const min_smap_sz(have_buildings() ? 1024U : 512U); // 1024x1024 seems to be required to prevent shadow artifacts on the sides of buildings
	unsigned lod_level(0);
	for (unsigned smap_sz = shadow_map_sz; (smap_sz > min_smap_sz && lod_level+1 < NUM_SMAP_LODS); ++lod_level) {smap_sz >>= 1;}
	return lod_level;
}
void tile_t::setup_shadow_maps(tile_shadow_map_manager &smap_manager, bool cleanup_only) {

	if (!shadow_map_enabled()) return; // disabled
	if (draw_model != 0) {clear_shadow_map(&smap_manager); return;} // skip shadow calculation in wireframe mode
	//timer_t timer("Create Tile Shadow Maps");
	float const smap_dist_scale(get_dist_to_camera_in_tiles(1)/(SMAP_NEW_THRESH*smap_thresh_scale));

	if (smap_dist_scale < 1.0) { // allocate new shadow maps or change shadow map LOD levels
		unsigned const max_lod_level(calc_max_smap_lod());
		float const lod_level_f(min(5.0f*smap_dist_scale, float(max_lod_level))); // clamp to max supported LOD
		unsigned lod_level(floor(lod_level_f)), lod_level_reduce(floor(max(0.0f, (lod_level_f - 0.1f)))); // add hysteresis
		
		// for very large shadow maps, limit LOD 0 to only the tile containing the camera
		if (shadow_map_sz >= 8192 && !get_mesh_bcube().contains_pt_xy(get_camera_pos())) {
			if (lod_level == 0) {lod_level_reduce = 1;} // no hysteresis for LOD 0
			lod_level = min(lod_level+1, max_lod_level);
		}
		if (lod_level_reduce > smap_lod_level) {clear_shadow_map(&smap_manager);} // LOD decrease
		if (cleanup_only) return; // done
		if (lod_level < smap_lod_level) {clear_shadow_map(&smap_manager);} // LOD increase

		if (smap_data.empty()) {
			for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) { // uses tu_id 13 and 14
				smap_data.push_back(smap_manager.new_smap_data(TILE_SMAP_START_TU_ID+i, this, i, lod_level));
			}
			//cout << "*** LOD: " << lod_level << " curLOD: " << smap_lod_level << " LODf: " << lod_level_f << " dist: " << smap_dist_scale << endl;
			smap_lod_level = lod_level;
		}
	}
	if (cleanup_only) return; // done
	cube_t bcube(get_shadow_bcube());
	// extend bcube upwards to include any models above the mesh that cast shadows on this tile
	// FIXME: still not correct for low sun pos - need a more accurate way to determine which models can shadow this tile
	cube_t const models_bcube(calc_and_return_all_models_bcube());
	if (models_bcube != all_zeros_cube) {max_eq(bcube.z2(), models_bcube.z2());}
	smap_data.create_if_needed(bcube);
}

bool tile_t::shadow_maps_allocated() const {
	for (unsigned i = 0; i < smap_data.size(); ++i) {if (smap_data[i].is_allocated()) return 1;}
	return 0;
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
	setup_texture(height_tid, 0, 0, 0, 0, 0);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, zvsize);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, stride, stride, 0, GL_RED, GL_FLOAT, &zvals.front());
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0); // reset to 0
}


void get_texture_ixs(int &sand_tex_ix, int &dirt_tex_ix, int &grass_tex_ix, int &rock_tex_ix, int &snow_tex_ix) {

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
}


bool check_region_int(cube_t const &region, vect_cube_t const &cubes) { // has_bcube_int_xy(), but without pad
	for (cube_t const &c : cubes) {
		if (c.intersects_xy(region)) return 1;
	}
	return 0;
}
void tile_t::create_texture(mesh_xy_grid_cache_t &height_gen) {

	//highres_timer_t timer("Create Tile Weights Texture"); // 1.38ms base, 1.5ms with buildings/roads/driveways/porches/doorsteps
	assert(zvals.size() == zvsize*zvsize);
	unsigned const tsize(stride), num_texels(tsize*tsize);
	int sand_tex_ix(-1), dirt_tex_ix(-1), grass_tex_ix(-1), rock_tex_ix(-1), snow_tex_ix(-1);
	get_texture_ixs(sand_tex_ix, dirt_tex_ix, grass_tex_ix, rock_tex_ix, snow_tex_ix);

	if (weight_tid == 0) { // create weights
		has_any_grass = has_tunnel = 0;
		grass_blocks.clear();
		mesh_weight_data.resize(4*num_texels); // RGBA
		unsigned const grass_block_dim(get_grass_block_dim());
		float const xy_mult(1.0/float(size)), water_level(get_water_z_height());
		float const MESH_NOISE_SCALE = 0.003;
		float const MESH_NOISE_FREQ  = 80.0;
		float const dz_inv(1.0f/(zmax - zmin));
		float const noise_scale(((mesh_gen_shape == 2) ? 2.0 : 1.0)*MESH_NOISE_SCALE*mesh_scale_z); // add more noise for ridged
		float const steep_mult_grass(1.0f/(sthresh[0][1] - sthresh[0][0]));
		float const steep_mult_snow (1.0f/(sthresh[1][1] - sthresh[1][0]));
		float const steep_mult_rock (1.0f/(0.8f*sthresh[0][0] - 0.5f*sthresh[0][0]));
		float const vnz_scale((mesh_gen_mode == MGEN_DWARP_GPU) ? SQRT2 : 1.0); // allow for steeper slopes when domain warping is used
		int const llc_x(x1 - xoff2), llc_y(y1 - yoff2);
		point const query_pos(get_xval(tsize/2 + llc_x), get_yval(tsize/2 + llc_y), 0.0); // in local tile space, not camera space
		bool const check_mesh_mask(check_mesh_disable(query_pos, radius)), check_buildings(no_grass_under_buildings());
		int k1, k2, k3, k4;
		height_gen.build_arrays(MESH_NOISE_FREQ*get_xval(x1), MESH_NOISE_FREQ*get_yval(y1), MESH_NOISE_FREQ*deltax,
			MESH_NOISE_FREQ*deltay, tsize, tsize, 0, 1); // force_sine_mode=1
		vector<float> rand_vals(tsize*tsize);
		bool row_ec_valid(0);
		vect_cube_t exclude_cubes, row_exclude_cubes, allow_cubes; // in camera space
		cube_t const mesh_bcube(get_mesh_bcube());
		get_city_grass_coll_cubes(mesh_bcube, exclude_cubes, allow_cubes);
		has_tunnel |= tile_contains_tunnel(mesh_bcube);

#pragma omp parallel for schedule(static,1) num_threads(2)
		for (int y = 0; y < (int)tsize-DEBUG_TILE_BOUNDS; ++y) {
			for (unsigned x = 0; x < tsize-DEBUG_TILE_BOUNDS; ++x) {
				rand_vals[y*tsize + x] = noise_scale*height_gen.eval_index(x, y, 50);
			}
		}
		for (unsigned y = 0; y < tsize-DEBUG_TILE_BOUNDS; ++y) { // not threadsafe
			float const yv(float(y)*xy_mult), ry(get_yval(y + llc_y + yoff)), radius_y(0.75*DY_VAL), ry1(ry - radius_y), ry2(ry + radius_y);
			row_ec_valid = 0;

			for (unsigned x = 0; x < tsize-DEBUG_TILE_BOUNDS; ++x) {
				unsigned const ix_val(y*tsize + x), off(4*ix_val);
				
				if (check_mesh_mask || inside_city == 1) { // have tunnels, or partially inside a city
					point const query_pos(get_xval(x + llc_x)+0.5*DX_VAL, get_yval(y + llc_y)+0.5*DY_VAL, 0.0); // global space
					
					if ((check_mesh_mask && check_mesh_disable(query_pos, HALF_DXY)) || (inside_city == 1 && check_inside_city(query_pos, HALF_DXY))) {
						mesh_weight_data[off+0] = mesh_weight_data[off+1] = 255; // set invalid values to flag as transparent
						mesh_weight_data[off+2] = mesh_weight_data[off+3] = 0;   // make sure grass is disabled
						has_tunnel = 1; // Note: should be covered by the tile_contains_tunnel(), but we include this case for safety
						continue;
					}
				}
				float weights[NTEX_DIRT] = {0};
				unsigned const ix(y*zvsize + x);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const mhmin(min(min(mh00, mh01), min(mh10, mh11))), mhmax(max(max(mh00, mh01), max(mh10, mh11)));
				float const rand_offset(rand_vals[y*tsize + x]);
				float const relh1(relh_adj_tex + (mhmin - zmin)*dz_inv + rand_offset), relh2(relh_adj_tex + (mhmax - zmin)*dz_inv + rand_offset);
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
					float vnz(vnz_scale*normal.z/normal.mag());
					// add random noise here as well to produce dry patches of dirt and sand in the grass
					if (grass && vnz > sti[1]) {vnz = CLIP_TO_01(1.0f + 20.0f*rand_offset);}

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
				float const xv(float(x)*xy_mult);

				// convert dirt to sand only when there is vegetation; even though it doesn't make sense to have dirt when there's no vegetation, it adds more texture variety
				if (vegetation > 0.0) {
					float const dirt_scale(BILINEAR_INTERP(params, dirt, xv, yv)); // slow

					if (dirt_scale < 1.0) { // apply dirt scale: convert dirt to sand
						weights[sand_tex_ix] += (1.0 - dirt_scale)*weights[dirt_tex_ix];
						weights[dirt_tex_ix] *= dirt_scale;
					}
				}
				if (grass) {
					float grass_scale((mhmin < water_level) ? 0.0f : BILINEAR_INTERP(params, grass, xv, yv)); // no grass under water
					bool replace_grass_with_dirt(0);

					if (grass_scale > 0.0 && !exclude_cubes.empty()) { // exclude bridges and tunnels here
						if (!row_ec_valid) { // optimization: calculate exclude cubes for the current yval row; not needed for allow_cubes because they're sparse
							cube_t const row_region(mesh_bcube.x1(), mesh_bcube.x2(), ry1, ry2, 0.0, 0.0);
							row_exclude_cubes.clear();

							for (cube_t const &c : exclude_cubes) {
								if (c.intersects_xy(row_region)) {row_exclude_cubes.push_back(c);}
							}
							row_ec_valid = 1;
						}
						if (!row_exclude_cubes.empty()) {
							float const rx(get_xval(x + llc_x + xoff)), radius_x(0.75*DX_VAL), rx1(rx - radius_x), rx2(rx + radius_x);
							cube_t const grass_region(rx1, rx2, ry1, ry2, 0.0, 0.0);
							replace_grass_with_dirt = (check_region_int(grass_region, row_exclude_cubes) && !check_region_int(grass_region, allow_cubes));
						}
					}
					if (!replace_grass_with_dirt && check_buildings && grass_scale > 0.0 && mh01 == mh00 && mh10 == mh00 && mh11 == mh00) { // look for area flattened under a building
						point const test_pt(get_xval(x + llc_x + xoff), ry, mh00); // in camera space
						replace_grass_with_dirt = check_buildings_no_grass(test_pt); // xy_only 1.61 => 1.76
					}
					if (replace_grass_with_dirt) {
						weights[dirt_tex_ix] += weights[grass_tex_ix]; // replace grass with dirt
						weights[grass_tex_ix] = 0.0;
						grass_scale = 0.0;
					}
					else if (grass_scale < 1.0) { // apply grass scale: convert grass to sand
						float const gscale(CLIP_TO_01(2.5f*(grass_scale - 0.5f) + 0.5f));
						weights[sand_tex_ix ] += (1.0 - gscale)*weights[grass_tex_ix];
						weights[grass_tex_ix] *= gscale;
					}
					if (grass_scale > 0.0) {add_grass_block_at(x, y, mhmin, mhmax, grass_block_dim);}
				} // end grass
				for (unsigned i = 0; i < NTEX_DIRT-1; ++i) { // Note: weights should sum to 1.0, so we can calculate w4 as 1.0-w0-w1-w2-w3
					mesh_weight_data[off+i] = ((weights[i] <= 0.01) ? 0 : ((weights[i] >= 0.99) ? 255 : (unsigned char)(255.0*weights[i])));
				}
			} // for x
		} // for y
	}
	else { // use existing weights
		assert(recalc_tree_grass_weights); // can only get here in this case
		assert(mesh_weight_data.size() == 4*num_texels);
	}
	weight_data = mesh_weight_data; // deep copy so that tree_map doesn't alter original weights

	if (!tree_map.empty()) {
		for (unsigned i = 0; i < num_texels; ++i) { // replace grass under trees with dirt
			unsigned const off(4*i);
			if (tree_map[i].ao == 255 || weight_data[off+grass_tex_ix] == 0) continue; // no trees or no grass
			float const v(tree_map[i].ao/255.0), w(weight_data[off+dirt_tex_ix] + (1.0 - v)*weight_data[off+grass_tex_ix]);
			weight_data[off+dirt_tex_ix]   = (unsigned char)(max(0.0f, min(255.0f, w)));
			weight_data[off+grass_tex_ix] *= v;
		}
	}
	recalc_tree_grass_weights = 0;
	create_or_update_weight_tex();
	calc_avg_mesh_color();
}


void tile_t::add_grass_block_at(unsigned x, unsigned y, float mhmin, float mhmax, unsigned grass_block_dim) {

	if (is_distant || x >= size || y >= size || !gen_grass_map()) return;
	unsigned const bx(x/GRASS_BLOCK_SZ), by(y/GRASS_BLOCK_SZ), bix(by*grass_block_dim + bx);
	if (grass_blocks.empty()) {grass_blocks.resize(grass_block_dim*grass_block_dim);}
	assert(bix < grass_blocks.size());
	grass_block_t &gb(grass_blocks[bix]);

	if (gb.ix == 0) { // not yet set
		gb.ix = (((x1 + x) + 1567*(y1 + y)) % num_rnd_grass_blocks) + 1; // select a random block
		//gb.ix = (rand() % num_rnd_grass_blocks) + 1; // select a random block
		gb.zmin = mhmin;
		gb.zmax = mhmax;
	}
	else {
		gb.zmin = min(gb.zmin, mhmin);
		gb.zmax = max(gb.zmax, mhmax);
	}
}


void tile_t::create_or_update_weight_tex() {

	unsigned const tsize(stride), num_texels(tsize*tsize);
	assert(weight_data.size() == 4*num_texels);

	if (weight_tid == 0) { // create weight texture
		setup_texture(weight_tid, 0, 0, 0, 0, 0);
		assert(weight_tid > 0 && glIsTexture(weight_tid));
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tsize, tsize, 0, GL_RGBA, GL_UNSIGNED_BYTE, &weight_data.front()); // internal_format = GL_COMPRESSED_RGBA - too slow
	}
	else { // update texture
		bind_2d_texture(weight_tid);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tsize, tsize, GL_RGBA, GL_UNSIGNED_BYTE, &weight_data.front());
	}
}


void tile_t::calc_avg_mesh_color() {

	// determine average mesh color - doesn't work well due to discontinuities at the tile boundaries, would be better if adjacent tile was known
	int tot_weights[NTEX_DIRT] = {0}; // use a signed value to avoid wrapping around to const max when there's an error
	avg_mesh_tex_color = colorRGB(0,0,0); // black
	unsigned const tsize(stride), num_texels(tsize*tsize);

	for (unsigned i = 0; i < num_texels; ++i) {
		unsigned const off(4*i);
		if (weight_data[off+1] == 255 && weight_data[off+1] == 255) continue; // skip disabled mesh sections
		for (unsigned j = 0; j < 4; ++j) {tot_weights[j] += weight_data[off+j];}
	}
	tot_weights[4] += (255*num_texels - tot_weights[0] - tot_weights[1] - tot_weights[2] - tot_weights[3]);
	for (unsigned i = 0; i < NTEX_DIRT; ++i) {avg_mesh_tex_color += get_avg_color_for_landscape_tex(i) * (tot_weights[i] / (255.0*num_texels));}
	//avg_mesh_tex_color *= (1.0/avg_mesh_tex_color.get_luminance());
	avg_mesh_tex_color *= 5.0;
	avg_mesh_tex_color  = (avg_mesh_tex_color*0.2 + WHITE*0.8); // 20% tinted
}


bool tile_t::update_range(tile_shadow_map_manager &smap_manager) { // if returns 0, tile will be deleted

	update_pine_tree_state(0); // can free pine tree vbos
	update_animals(); // if any were generated
	float const dist(get_rel_dist_to_camera());
	
	if (dist > CLEAR_DIST_TILES || mesh_height_invalid) {
		if (!just_cleared) {clear_vbo_tid(&smap_manager);} // avoid clearing every frame
		just_cleared = 1;
	}
	else {just_cleared = 0;}
	if (dist*TILE_RADIUS > SMAP_DEL_THRESH*smap_thresh_scale) {clear_shadow_map(&smap_manager);} // too far, delete old shadow maps
	return (dist < DELETE_DIST_TILES && !mesh_height_invalid);
}


// *** trees ***

bool tile_t::can_have_pine_palm_trees() const {
	return (pine_trees_enabled() && can_have_trees() && can_have_pine_palm_trees_in_zrange(mzmin, mzmax, 1)); // skip_range_check_if_manually_placed=1
}
bool tile_t::can_have_decid_trees() const {
	return (decid_trees_enabled() && can_have_trees() && can_have_decid_trees_in_zrange   (mzmin, mzmax, 1)); // skip_range_check_if_manually_placed=1
}

void tile_t::init_pine_tree_draw() {

	float const density[4] = {params[0][0].veg, params[0][1].veg, params[1][0].veg, params[1][1].veg};
	ptree_off.set_from_xyoff2();
	if (enable_instanced_pine_trees()) {pine_trees.enable_instanced();}
	pine_trees.gen_trees(x1+ptree_off.dxoff, y1+ptree_off.dyoff, x2+ptree_off.dxoff, y2+ptree_off.dyoff, density);
	postproc_trees(pine_trees, ptzmax);
}

void tile_t::update_pine_tree_state(bool upload_if_needed, bool force_high_detail) {

	if (!pine_trees_enabled() || pine_trees.empty()) return;
	float const weight(get_tree_far_weight(force_high_detail, (pine_trees.num_palm_trees > 0))); // force high detail trees during the shadow map pass
	float const weights[2] = {1.0f-weight, weight}; // {high, low} detail
	//timer_t timer("Update Pine Trees");

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
void tile_t::draw_tree_leaves_lod(shader_t &s, vector3d const &xlate, bool low_detail, int xlate_loc) {
	bool const draw_all(low_detail || camera_pdu.point_visible_test(get_center())); // tile center is in view
	pine_trees.draw_pine_leaves(s, 0, low_detail, draw_all, contains_camera(), xlate, xlate_loc);
}

// Note: xlate has different meanings here: for near leaves, it's a vec3 attribute; for far leaves, it's a vec2 uniform; for branches, it's -1
void tile_t::draw_pine_trees(shader_t &s, vector<tile_t *> &to_draw_trunk_pts, bool draw_trunks, bool draw_near_leaves,
	bool draw_far_leaves, bool shadow_pass, bool reflection_pass, bool enable_smap, int xlate_loc)
{
	if (pine_trees.empty() || /*!can_have_pine_palm_trees()*/!can_have_trees()) return; // Note: skip water check for palm trees
	//timer_t timer("Draw Pine Trees");
	point const camera(get_camera_pos());
	vector3d const xlate(ptree_off.get_xlate());
	vector2d const camera_xlate(xlate.x-camera.x, xlate.y-camera.y);
	bool const force_high_detail(shadow_pass);
	fgPushMatrix();
	translate_to(xlate);
	
	if (draw_trunks) {
		float const dscale(shadow_pass ? 0.0 : get_tree_dist_scale()); // force polygons in shadow pass

		if (dscale < 1.0) { // close, draw as polygons
			if (enable_smap) {bind_and_setup_shadow_map(s); enable_smap = 0;}
			if (!shadow_pass) {set_mesh_ambient_color(s);}
			bool const all_visible(camera_pdu.sphere_visible_test(get_center(), -0.5*radius));
			if (!pine_trees.draw_trunks(shadow_pass, all_visible, 1, xlate)) {to_draw_trunk_pts.push_back(this);} // skip_lines=1; if some lines were skipped, add trunk pts to draw
		}
		else if (dscale < 2.0 && get_tree_far_weight(force_high_detail) < 0.5) { // far away, use low detail branches
			to_draw_trunk_pts.push_back(this);
		} // else very far, skip branches
	}
	if (draw_near_leaves || draw_far_leaves) { // could use reflection_pass as an optimization
		bool const has_palm(pine_trees.num_palm_trees > 0);
		float const weight(1.0 - get_tree_far_weight(force_high_detail, has_palm)); // 0 => low detail, 1 => high detail
		
		if (draw_far_leaves && weight < 1.0) {
			assert(xlate_loc >= 0);
			s.set_uniform_vector2d(xlate_loc, camera_xlate);
		}
		if (weight > 0.0 && weight < 1.0) { // use geomorphing with dithering (since alpha doesn't blend in the correct order)
			if (draw_near_leaves) {
				if (enable_smap) {bind_and_setup_shadow_map(s); enable_smap = 0;}
				int const loc(s.get_uniform_loc("max_noise"));
				s.set_uniform_float(loc, weight);
				draw_tree_leaves_lod(s, xlate, 0, xlate_loc); // near leaves
				s.set_uniform_float(loc, 1.0);
			}
			if (draw_far_leaves) {
				int const loc(s.get_uniform_loc("min_noise"));
				s.set_uniform_float(loc, weight);
				draw_tree_leaves_lod(s, xlate, 1, xlate_loc); // far leaves
				s.set_uniform_float(loc, 0.0);
			}
		}
		else if ((weight == 0.0) ? draw_far_leaves : draw_near_leaves) {
			if (enable_smap) {bind_and_setup_shadow_map(s); enable_smap = 0;}
			draw_tree_leaves_lod(s, xlate, (weight == 0.0), xlate_loc);
		}
		if (draw_near_leaves && has_palm && (weight > 0.0 || shadow_pass || (pine_trees.instanced && get_tree_dist_scale(has_palm) < 3.0))) { // draw palm trees
			if (enable_smap) {bind_and_setup_shadow_map(s); enable_smap = 0;}
			s.add_uniform_int("use_bent_quad_tcs", 1);
			if (pine_trees.instanced) {pine_trees.draw_palm_insts(s, camera_pdu.sphere_visible_test(get_center(), -0.5*radius), xlate, xlate_loc);}
			else {pine_trees.draw_non_pine_leaves(0, 1, 0, xlate_loc, -1, xlate);}
			s.add_uniform_int("use_bent_quad_tcs", 0);
		}
	}
	fgPopMatrix();
}

void tile_t::draw_trunk_pts(shader_t &s) {

	assert(!pine_trees.empty());
	s.add_uniform_vector3d("xlate", ptree_off.get_xlate());
	pine_trees.draw_trunk_pts();
}


void tile_t::gen_decid_trees_if_needed() {

	if (num_trees > 0 && max_unique_trees == 0) {
		cout << "Warning: max_unique_trees needs to be set to something reasonable for tiled terrain mode trees to work efficiently. Setting to 100." << endl;
		max_unique_trees = 100;
	}
	if (decid_trees.was_generated() || !can_have_decid_trees()) return; // already generated, elevation too high, or distant tile (no trees yet)
	assert(decid_trees.empty());
	dtree_off.set_from_xyoff2();
	decid_trees.gen_deterministic(x1+dtree_off.dxoff, y1+dtree_off.dyoff, x2+dtree_off.dxoff, y2+dtree_off.dyoff, vegetation*get_avg_veg(), mesh_dz, this);
	postproc_trees(decid_trees, dtzmax);
}

void tile_t::update_decid_trees() {
	if (decid_trees_enabled()) {decid_trees.check_render_textures();}
}

void tile_t::set_mesh_ambient_color(shader_t &s) const {s.add_uniform_color("ambient_tint", avg_mesh_tex_color);}

void tile_t::draw_decid_trees(shader_t &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves, bool reflection_pass, bool shadow_pass, bool enable_smap) {

	if (decid_trees.empty() || !can_have_trees()) return;
	//timer_t timer("Draw Decid Trees");
	if (enable_smap) {bind_and_setup_shadow_map(s);}
	if (draw_branches && !shadow_pass) {set_mesh_ambient_color(s);}
	// Note: shadow_only mode doesn't help performance much
	decid_trees.draw_branches_and_leaves(s, lod_renderer, draw_branches, draw_leaves, shadow_pass, reflection_pass, dtree_off.get_xlate());
}


// *** scenery/grass/flowers ***

void tile_t::update_scenery() {

	if (!scenery_enabled() || is_distant) return; // no scenery
	float const dist_scale(get_scenery_dist_scale(0)); // tree_dist_scale should correlate with mesh scale
	if (scenery.generated && dist_scale > 1.2) {scenery.clear();} // too far away
	if (scenery.generated || dist_scale > 1.0 || !is_visible()) return; // already generated, too far away, or not visible
	//timer_t timer("Gen Scenery");
	scenery_off.set_from_xyoff2();
	// set our tile's decid trees as current so that logs and stumps get the correct colors; assumes trees are generated before scenery
	scenery.gen(x1+scenery_off.dxoff, y1+scenery_off.dyoff, x2+scenery_off.dxoff, y2+scenery_off.dyoff, vegetation*get_avg_veg(), 1, decid_trees);
}

void tile_t::draw_scenery(shader_t &s, shader_t &vrs, bool draw_opaque, bool draw_leaves, bool reflection_pass, bool shadow_pass, bool enable_shadow_maps) {

	if (!scenery.generated || get_scenery_dist_scale(reflection_pass) > 1.0) return;
	//timer_t timer("Draw Scenery");
	fgPushMatrix();
	vector3d const xlate(scenery_off.get_xlate());
	translate_to(xlate);
	float const scale_val(get_scenery_thresh(reflection_pass)*get_tile_width());
	if (enable_shadow_maps) {bind_and_setup_shadow_map(s);}
	if (draw_opaque) {scenery.draw_opaque_objects(s, vrs, shadow_pass, xlate, 0, scale_val, reflection_pass);} // shader not passed in here
	if (draw_leaves) {scenery.draw_plant_leaves  (s, shadow_pass, xlate, reflection_pass);}
	fgPopMatrix();
}

void tile_t::pre_draw_grass_flowers(shader_t &s, bool use_cloud_shadows) const {

	bind_texture_tu(height_tid, 2);
	bind_texture_tu(normal_tid, 4);
	bind_texture_tu(shadow_tid, 6);
	bind_and_setup_shadow_map(s);
	if (use_cloud_shadows) {s.add_uniform_vector3d("cloud_offset", vector3d(get_xval(x1), get_yval(y1), 0.0));}
	s.add_uniform_vector2d("xlate", vector2d(get_xval(x1 + xoff - xoff2), get_yval(y1 + yoff - yoff2)));
}

unsigned tile_t::draw_grass(shader_t &s, vector<vector<vector2d> > *insts, bool use_cloud_shadows, bool enable_tess, int lt_loc) {

	if (!has_grass()) return 0; // or can test has_any_grass
	float const grass_thresh(get_grass_thresh_pad());
	point const camera(get_camera_pos());
	if (get_min_dist_to_pt(camera) > grass_thresh) return 0; // too far away to draw
	pre_draw_grass_flowers(s, use_cloud_shadows);
	bind_texture_tu(weight_tid, 3);
	unsigned const grass_block_dim(get_grass_block_dim());
	assert(grass_blocks.size() == grass_block_dim*grass_block_dim);
	float const llcx(get_xval(x1 + xoff - xoff2)), llcy(get_yval(y1 + yoff - yoff2));
	float const dx_step(GRASS_BLOCK_SZ*deltax), dy_step(GRASS_BLOCK_SZ*deltay);
	float const lod_scale(GRASS_LOD_SCALE/(tt_grass_scale_factor*get_scaled_tile_radius()));
	float const block_grass_thresh(grass_thresh + (SQRT2*radius)/grass_block_dim), bg_thresh_sq(block_grass_thresh*block_grass_thresh);
	point const adj_camera(camera + point(0.0, 0.0, 2.0*grass_length));
	bool const all_visible(camera_pdu.cube_completely_visible(get_mesh_bcube()));
	unsigned num_drawn(0);

	for (unsigned y = 0; y < grass_block_dim; ++y) {
		for (unsigned x = 0; x < grass_block_dim; ++x) {
			grass_block_t const &gb(grass_blocks[y*grass_block_dim+x]);
			if (gb.ix == 0) continue; // empty block
			float const bcx1(llcx + x*dx_step), bcy1(llcy + y*dy_step);
			cube_t const bcube(bcx1, bcx1+dx_step, bcy1, bcy1+dy_step, gb.zmin, (gb.zmax + grass_length));
			float const dist_sq(p2p_dist_sq(camera, bcube.closest_pt(camera)));
			if (dist_sq > bg_thresh_sq || (!all_visible && !camera_pdu.cube_visible(bcube))) continue;

			if (dist_sq < 0.56*bg_thresh_sq) { // only do back face culling on nearby blocks (75% of max dist)
				bool back_facing(1);

				for (unsigned yy = y*GRASS_BLOCK_SZ; yy <= (y+1)*GRASS_BLOCK_SZ && back_facing; ++yy) { // Note: could precompute avg normal per block, but it probably isn't necessary
					for (unsigned xx = x*GRASS_BLOCK_SZ; xx <= (x+1)*GRASS_BLOCK_SZ && back_facing; ++xx) {
						unsigned const ix(yy*zvsize + xx);
						back_facing &= (dot_product(get_norm_not_normalized(ix), (adj_camera - point(llcx+xx*deltax, llcy+yy*deltay, zvals[ix]))) < 0.0);
					}
				}
				if (back_facing) continue;
			}
			unsigned const lod_level(min(NUM_GRASS_LODS-1, unsigned(lod_scale*sqrt(dist_sq))));
			assert(insts != NULL);
			insts[lod_level].resize(num_rnd_grass_blocks); // may already be the correct size
			unsigned const bix(gb.ix - 1);
			assert(bix < num_rnd_grass_blocks);
			insts[lod_level][bix].emplace_back(x*dx_step, y*dy_step);
		} // for x
	} // for y
	for (unsigned lod = 0; lod < NUM_GRASS_LODS; ++lod) {
		for (unsigned bix = 0; bix < insts[lod].size(); ++bix) {
			vector<vector2d> &v(insts[lod][bix]);
			if (v.empty()) continue;
			glVertexAttribPointer(lt_loc, 2, GL_FLOAT, GL_FALSE, sizeof(vector2d), get_dynamic_vbo_ptr(&v.front(), v.size()*sizeof(vector2d)));
			num_drawn += grass_tile_manager.render_block(bix, lod, 1.0, v.size(), enable_tess);
			v.clear();
		} // for bix
	} // for lod
	return num_drawn;
}

unsigned tile_t::draw_flowers(shader_t &s, bool use_cloud_shadows) {

	if (!has_grass()) return 0; // no grass, no flowers
	float const flower_thresh(FLOWER_REL_DIST*get_grass_thresh_pad());
	if (get_min_dist_to_pt(get_camera_pos()) > flower_thresh) return 0; // too far away to draw
	flowers.gen_flowers(weight_data, stride, x1-xoff2, y1-yoff2); // mesh weight + tree dirt
	if (flowers.empty()) return 0; // no flowers generated
	pre_draw_grass_flowers(s, use_cloud_shadows);
	flowers.check_vbo();
	flowers.draw_triangles(s);
	return flowers.size();
}

bool tile_t::choose_butterfly_dest(point &dest, sphere_t &plant_bsphere, rand_gen_t &rgen) const {
	if (!scenery.generated) return 0; // no scenery
	if (!scenery.choose_butterfly_dest(dest, plant_bsphere, rgen)) return 0;
	vector3d const xlate(scenery_off.get_xlate() - get_camera_coord_space_xlate());
	dest += xlate; // convert to world space
	plant_bsphere.pos += xlate;
	return 1;
}


// *** clouds ***

void tile_cloud_t::draw(vpc_shader_t &s, vector3d const &xlate, float alpha_mult) const {

	vector3d const view_dir(get_camera_pos() - (pos + xlate));
	float const view_dist(view_dir.mag()), max_dist(max(get_draw_tile_dist(), get_inf_terrain_fog_dist()));
	if (view_dist >= max_dist) return; // too distant to draw
	float const val(1.0 - (max_dist - view_dist)/max_dist), alpha(alpha_mult*(1.0f - val*val));
	colorRGBA const cloud_color(get_cloud_color());
	s.set_uniform_color(s.c1i_loc, colorRGBA(cloud_color*0.75, alpha)); // inner color
	s.set_uniform_color(s.c1o_loc, colorRGBA(cloud_color*1.00, alpha)); // outer color
	s.set_uniform_float(s.rad_loc, get_rmax());
	s.set_uniform_float(s.off_loc, (pos_hash + 0.002*tfticks)); // used as a hash
	s.set_uniform_vector3d(s.vd_loc, view_dir/view_dist); // local object space
	s.set_uniform_vector3d(s.rs_loc, size/get_rmax());
	fgPushMatrix();
	translate_to(pos + xlate);
	draw_quads(1); // depth map is disabled in the caller
	fgPopMatrix();
}

void tile_cloud_manager_t::gen_new_cloud() {

	push_back(tile_cloud_t());
	tile_cloud_t &c(back());
	c.pos   = rgen.gen_rand_cube_point(range);
	c.size  = vector3d(rgen.rand_uniform(1.0, 2.0), rgen.rand_uniform(1.0, 2.0), rgen.rand_uniform(0.6, 1.0)); // smaller in z
	c.size *= rgen.rand_uniform(2.0, 4.0);
	c.pos_hash = c.pos.x;
	c.gen_pts(c.size);
	update_bcube(c);
}

void tile_cloud_manager_t::gen(int x1, int y1, int x2, int y2) {

	if (generated) return; // already generated
	generated = 1;
	rgen.set_state(x1+123, y1+321);
	float const z_range(zmax - zmin), z_cloud_bot(zmax + cloud_height_offset*z_range + 2.0); // shift up slightly so that bottom of cloud is here
	range = cube_t(get_xval(x1), get_xval(x2), get_yval(y1), get_yval(y2), (z_cloud_bot + 0.0*z_range), (z_cloud_bot + 0.9*z_range));
	choose_num_clouds();
	populate_clouds();
}

void tile_cloud_manager_t::choose_num_clouds() { // this tile will always create a constant number of clouds
	num_clouds = ((clouds_per_tile == 0.0) ? 0 : unsigned(max(0.0f, rgen.rand_gaussian(clouds_per_tile, 4.0))));
}

void tile_cloud_manager_t::populate_clouds() {

	assert(empty());
	bcube.set_to_zeros();
	reserve(num_clouds);
	for (unsigned n = 0; n < num_clouds; ++n) {gen_new_cloud();}
}

void tile_cloud_manager_t::try_add_cloud(tile_cloud_t const &cloud) {

	if (!generated) return;
	//assert(range.contains_pt(cloud.pos)); // TESTING
	push_back(cloud);
	update_bcube(cloud);
}

void tile_cloud_manager_t::update_bcube(tile_cloud_t const &c) {
	if (bcube.is_all_zeros()) {bcube = c.get_bcube();}
	else {bcube.union_with_cube(c.get_bcube());}
}

void tile_cloud_manager_t::move_by_wind(tile_t const &tile) {

	if (!generated || wind == zero_vector) return; // not yet generated, or no wind
	vector3d const move_amt((0.01*fticks)*vector3d(wind.x, wind.y, 0.0)); // no z movement
	float const move_amt_mag(move_amt.mag());
	cur_move_dist += move_amt_mag;

	if (cur_move_dist > 2.0*get_tile_width()) { // when wind has moved the clouds an entire tile's width, choose a new num_clouds for edge tiles
		choose_num_clouds();
		cur_move_dist = 0.0;
	}
	if (empty() && num_clouds > 0) { // no current clouds
		// see if we're on the edge of the mesh that should generate incoming clouds
		bool on_edge(0);

		for (unsigned d = 0; d < 2; ++d) {
			int dxy[2] = {0,0};
			if (wind[d] < 0.0) {dxy[d] = 1;} else if (wind[d] > 0.0) {dxy[d] = -1;}
			if (tile.get_adj_tile(dxy[0], dxy[1]) == nullptr) {on_edge = 1; break;}
		}
		if (on_edge) {populate_clouds();} // create more random clouds
		return;
	}
	if (empty()) return;
	float const ws_mult(0.5/range.dz());
	bcube.set_to_zeros();

	for (auto i = begin(); i != end(); ++i) {
		float const wscale(1.5f - ws_mult*(i->pos.z - range.d[2][0])); // lower clouds have higher wind speed
		i->pos      += wscale*move_amt;
		i->pos_hash += 0.15*move_amt_mag;
		int dx(0), dy(0);
		if (i->pos.x < range.d[0][0]) {dx = -1;} else if (i->pos.x > range.d[0][1]) {dx = 1;}
		if (i->pos.y < range.d[1][0]) {dy = -1;} else if (i->pos.y > range.d[1][1]) {dy = 1;}
		if (dx == 0 && dy == 0) {update_bcube(*i); continue;} // still within the tile's bounds, keep it
		tile_t *adj(tile.get_adj_tile(dx, dy));
		if (adj != nullptr) {adj->add_cloud(*i);} // move to neighbor tile
		std::swap(*i, back()); pop_back(); --i; // remove this cloud
	} // for i
}

void tile_cloud_manager_t::get_draw_list(cloud_draw_list_t &clouds_to_draw, float mesh_zmin, float mesh_zmax) const {

	if (empty()) return;
	vector3d const xlate(get_tiled_terrain_model_xlate());

	for (auto i = begin(); i != end(); ++i) {
		point const pt(i->pos + xlate);
		if (!camera_pdu.sphere_visible_test(pt, i->get_rmax())) continue; // VFC
		float const zbot(i->pos.z - i->size.z);
		if (zbot < mesh_zmin) continue; // definitely below the mesh
		if (zbot + 0.2*i->size.z < water_plane_z) continue; // below the water
		float const start_dist(1.0*i->size.z);
		float alpha(1.0);

		if (zbot < mesh_zmax+start_dist) { // may be below the mesh
			float const dist(zbot - get_exact_zval(i->pos.x-xoff2*DX_VAL, i->pos.y-yoff2*DY_VAL));
			if (dist < 0.0) continue; // below the mesh
			if (dist < start_dist) {alpha = dist/start_dist;}
		}
		clouds_to_draw.push_back(cloud_inst_t(*i, -p2p_dist_sq(get_camera_pos(), pt), alpha));
	} // for i
}

unsigned tile_t::update_tile_clouds() {

	if (animate2) {clouds.move_by_wind(*this);}
	clouds.gen(x1, y1, x2, y2);
	return clouds.size();
}


// *** animals ***

template<typename A> void tile_t::propagate_animals_to_neighbor_tiles(animal_group_t<A> &animals) {

	if (animals.empty()) return;
	cube_t range(get_mesh_bcube_global());

	for (unsigned i = 0; i < animals.size(); ++i) {
		if (!animals[i].is_enabled()) {animals.remove(i); --i; continue;} // remove disabled animals
		point const &pos(animals[i].pos);
		int dx(0), dy(0);
		if      (pos.x < range.d[0][0]) {dx = -1;} // move -x
		else if (pos.x > range.d[0][1]) {dx =  1;} // move +x
		if      (pos.y < range.d[1][0]) {dy = -1;} // move -y
		else if (pos.y > range.d[1][1]) {dy =  1;} // move +y
		if (dx == 0 && dy == 0) continue; // position is within the current tile, keep this animal here
		tile_t *adj(get_adj_tile(dx, dy));
		if (adj != NULL) {adj->add_animal(animals[i]);} // move to adjacent tile, if one is present - pos should be valid within that tile
		animals.remove(i); --i; // remove from current tile
	}
}

void tile_t::update_animals() {

	if (!ENABLE_ANIMALS) return;
	// Note: animals are in global space, rather than camera local space like everything else;
	// this means that there will be FP errors when the player is far from the origin - but since fish and birds are "distant" objects, that may be okay;
	// using camera space is more difficult due to all of update code interacting with the rest of the scene, which is in global space;
	// also, animals can move between tiles, which complicates the math (since adjacent tiles may create their animals in a different starting camera space)
	if (water_is_lava || temperature > 0.5*WATER_MAX_TEMP || vegetation < 0.1) {} // no fish
	else if (!fish.was_generated()) {
		cube_t range(get_mesh_bcube_global());
		range.z2() = water_plane_z; // z extends from lowest mesh point to water surface
		fish.gen(num_fish_per_tile, range, this); // Note: could use get_water_bcube() for tighter range
	}
	else {
		fish.update(this);
		propagate_animals_to_neighbor_tiles(fish);
	}
	if (atmosphere > 0.4 && vegetation > 0.2) { // need atmosphere and vegetation for birds and butterflies
		adj_tiles_t adj_tiles;

		if (!birds.was_generated()) {
			cube_t range(get_mesh_bcube_global());
			float const z_range(zmax - zmin);
			range.z1() = zmax;
			range.z2() = zmax + 0.50*z_range; // Note: may be in the clouds
			birds.gen(num_birds_per_tile, range, this);
		}
		else {
			birds.flock(this, adj_tiles);
			birds.update(this);
			propagate_animals_to_neighbor_tiles(birds);
		}
		if (weight_tid == 0) {} // weight texture not yet generated, don't add until we know if there's any grass
		else if (!bflies.was_generated()) {
			bflies.gen(num_bflies_per_tile, get_mesh_bcube_global(), this);
		}
		else {
			bflies.run_mating(this, adj_tiles);
			bflies.update(this);
			propagate_animals_to_neighbor_tiles(bflies);
		}
	}
}


// *** rendering ***


void tile_t::pre_draw(mesh_xy_grid_cache_t &height_gen) {

	assert(!zvals.empty());
	if (tree_map.empty() && any_trees_enabled()) {apply_tree_ao_shadows();}
	if (weight_tid == 0 || recalc_tree_grass_weights) {create_texture(height_gen);}
	check_shadow_map_and_normal_texture();
	ensure_height_tid();
}


unsigned tile_t::get_lod_level(bool reflection_pass) const {

	if (mzmax == mzmin) {return NUM_LODS-1;} // flat, use lowest LOD (city plots)
	unsigned lod_level((reflection_pass && min_normal_z > 0.1) ? min(NUM_LODS-1, 1U) : 0);
	float dist(get_dist_to_camera_in_tiles(0));
	if (shadow_maps_allocated()) {dist /= 2;} // decrease distance for higher LOD when shadow maps are in use to avoid shadow artifacts on low res mesh
	
	if (min_normal_z > 0.0) { // normals have been calculated, adjust detail based on max slope
		if (min_normal_z > 0.9) { // flat, lower detail
			if (!is_water_enabled() || mzmax < water_plane_z || mzmin > water_plane_z) { // use get_water_z_height()?
				dist *= (min_normal_z > 0.95) ? 4.0 : 2.0;
			}
		}
		else if (min_normal_z < 0.25) {dist /= -log2(2.0*min_normal_z);} // high slope, higher detail
	}
    while (dist > (reflection_pass ? 1.8 : 2.0) && lod_level+1 < NUM_LODS && (size>>(lod_level+1)) >= 4) { // never divide smaller than 4x4 quads
        dist /= 2;
        ++lod_level;
    }
	assert(lod_level < NUM_LODS);
	return lod_level;
}

void disable_shadow_maps(shader_t &s) { // Note: uses different TUs compared to bind_default_sun_moon_smap_textures()
	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid_and_enabled(l)) continue;
		bind_texture_tu(get_empty_smap_tid(), TILE_SMAP_START_TU_ID+l); // bind empty shadow map
		s.add_uniform_float((l ? "sm_scale1" : "sm_scale0"), -1.0); // set shadow map scale to 0
	}
}

void tile_t::shader_shadow_map_setup(shader_t &s) const {
	// Note: some part of this call is shared across all tiles; however, in the case where more than one smap light is enabled,
	// the tu_id and enables may alternate between values for each tile, requiring every uniform to be reset per tile anyway
	if (smap_data.empty()) { // disable shadow map lookup when shadow map textures are unavailable
		if (!shadow_maps_disabled) {disable_shadow_maps(s);} // optimization: skip if already disabled (using a global variable)
		shadow_maps_disabled = 1;
	}
	else {
		smap_data.set_for_all_lights(s, nullptr); // no MVM
		shadow_maps_disabled = 0;
	}
}
void tile_t::bind_and_setup_shadow_map(shader_t &s) const {
	if (shadow_map_enabled()) {shader_shadow_map_setup(s);}
}
bool tile_t::try_bind_shadow_map(shader_t &s, bool check_only, unsigned *lod_level) const {
	if (!shadow_map_enabled() || smap_data.empty()) return 0;
	if (get_dist_to_camera_in_tiles(1) > SMAP_FADE_THRESH*smap_thresh_scale) return 0; // too far to need smap, even if it exists
	if (!check_only) {smap_data.set_for_all_lights(s, nullptr);}
	
	if (lod_level) { // return lod_level of first enabled shadow map
		for (unsigned i = 0; i < smap_data.size(); ++i) {
			if (is_light_enabled(i)) {*lod_level = smap_data[i].lod_level; break;}
		}
	}
	return 1;
}


void tile_t::bind_textures() const {

	assert(weight_tid > 0);
	assert(height_tid > 0);
	assert(normal_tid > 0);
	assert(shadow_tid > 0);
	bind_2d_texture(weight_tid);
	bind_texture_tu(height_tid, 12);
	bind_texture_tu(normal_tid, 7);
	bind_texture_tu(shadow_tid, 15);
}


void crack_ibuf_t::gen_offsets(vector<unsigned> &indices, unsigned size) {

	offsets.clear();
	offsets.resize(4*NUM_LODS*NUM_LODS);

	for (unsigned dim = 0; dim < 2; ++dim) { // x,y
		for (unsigned dir = 0; dir < 2; ++dir) { // lo, hi
			for (unsigned cur_lod = 0; cur_lod < NUM_LODS; ++cur_lod) {
				for (unsigned adj_lod = 0; adj_lod < NUM_LODS; ++adj_lod) {
					ix_sz_pair &ixsz(offsets[get_index(dim, dir, cur_lod, adj_lod)]);
					ixsz.ix = (unsigned short)indices.size();
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
					ixsz.sz = (unsigned short)indices.size() - ixsz.ix; // often 0
				} // for adj_lod
			} // for cur_lod
		} // for dir
	} // for dim
}


unsigned crack_ibuf_t::get_index(unsigned dim, unsigned dir, unsigned cur_lod, unsigned adj_lod) const {
	assert(dim < 2 && dir < 2 && cur_lod < NUM_LODS && adj_lod < NUM_LODS);
	return (NUM_LODS*(NUM_LODS*(2*dim + dir) + cur_lod) + adj_lod);
}


void tile_t::draw_mesh_vbo(indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], unsigned lod_level) const {

	assert(size > 0);
	unsigned const num_ixs(ivbo_ixs[lod_level+1] - ivbo_ixs[lod_level]);
	vbo_mgr.pre_render();
	vert_wrap_t::set_vbo_arrays(0); // normals are stored in normal_tid, tex coords come from texgen, color is constant
	glDrawRangeElements(GL_TRIANGLE_STRIP, 0, stride*stride, num_ixs, GL_UNSIGNED_INT, (void *)(ivbo_ixs[lod_level]*sizeof(unsigned)));
	++num_frame_draw_calls;
}

void tile_t::draw(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], crack_ibuf_t const &crack_ibuf, int reflection_pass, int shader_locs[2]) const {

	//timer_t timer("Draw Tile Mesh");
	// check if the tile was visible in the building mirror reflection but not in normal view (so wasn't setup)
	//if (get_checkerboard_bit()) return; // checkerboard drawing, for debugging
	if (inside_city == 2) return; // don't need to draw terrain if fully inside a city, but still need to draw trees, etc.
	if (!(weight_tid > 0 && height_tid > 0 && normal_tid > 0 && shadow_tid > 0)) return; // textures not yet created
	fgPushMatrix();
	vector3d const xlate(get_mesh_xlate());
	translate_to(xlate); // Note: not easy to replace with a uniform, due to texgen and fog dist calculations in the shader
	bool const draw_near_water(!is_distant && !reflection_pass && is_water_enabled() && has_water());

	if (!reflection_pass && cloud_shadows_enabled()) {
		s.ensure_uniform_loc(shader_locs[0], "cloud_offset");
		s.set_uniform_vector3d(shader_locs[0], vector3d(get_xval(x1), get_yval(y1), 0.0));
	}
	if (draw_near_water) { // for underwater caustics texture
		s.ensure_uniform_loc(shader_locs[1], "tc_xlate");
		s.set_uniform_vector3d(shader_locs[1], xlate);
	}
	if (reflection_pass == 2) {disable_shadow_maps(s);} // disabled for mirror reflections because shadows don't work
	else {shader_shadow_map_setup(s);}
	bind_textures(); // Note: moved after the disable_shadow_maps() call to ensure TU 0 is not overwritten
	unsigned const lod_level(get_lod_level(reflection_pass));
	draw_mesh_vbo(vbo_mgr, ivbo_ixs, lod_level);
	
	// draw cracks
	for (unsigned dim = 0; dim < 2; ++dim) { // x,y
		for (unsigned dir = 0; dir < 2; ++dir) { // lo, hi
			int dx(0), dy(0);
			(dim ? dy : dx) += (dir ? 1 : -1);
			tile_t *adj(get_adj_tile(dx, dy)); // need to handle (is_distant != adj->is_distant) when distant tiles are implemented
			if (adj == NULL) continue; // no adjacent tile
			//if (!adj->is_visible() || !adj->rel_dist_to_camera_xy_lt(DRAW_DIST_TILES)) continue;
			ix_sz_pair const &ixsz(crack_ibuf.lookup(crack_ibuf.get_index(dim, dir, lod_level, adj->get_lod_level(reflection_pass))));
			if (ixsz.sz == 0) continue;
			glDrawRangeElements(GL_TRIANGLES, 0, stride*stride, ixsz.sz, GL_UNSIGNED_INT, (void *)(ixsz.ix*sizeof(unsigned)));
			++num_frame_draw_calls;
		} // for dir
	} // for dim
	bind_vbo(0, 0); // unbind vertex buffer
	fgPopMatrix();
	if (draw_near_water) {draw_water_cap(s, 1);}
}

void tile_t::draw_shadow_pass(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], int xlate_loc) { // not const because creates height_tid

	if (inside_city == 2) return; // don't need to draw terrain if fully inside a city
	ensure_height_tid();
	s.set_uniform_vector3d(xlate_loc, get_mesh_xlate());
	assert(height_tid > 0);
	bind_texture_tu(height_tid, 12);
	draw_mesh_vbo(vbo_mgr, ivbo_ixs, 0); // LOD is always 0
	bind_vbo(0, 0); // unbind vertex buffer
}


void tile_t::draw_water_cap(shader_t &s, bool textures_already_set) const {

	bool const dist_water(draw_distant_water());
	//if (dist_water) return; // skip water cap
	cube_t const bcube(get_mesh_bcube());
	static vector<vert_wrap_t> wverts; // reused across tiles/frames
	wverts.clear();

	// draw vertical edges that cap the water volume and will be blended between underwater black and fog colors
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			int dxy[2] = {0,0};
			dxy[dim] = (dir ? 1 : -1);
			tile_t const *const adj_tile(get_adj_tile(dxy[0], dxy[1]));
				
			if (adj_tile && !adj_tile->is_distant) {
				if (adj_tile->rel_dist_to_camera_xy_lt(DRAW_DIST_TILES)) continue;
				if (!adj_tile->has_water ()) continue; // adj tile has no water, so we can't have any uncapped water on this edge
				if (!adj_tile->is_visible()) continue; // adj tile not visible,  so we can't see this edge
			}
			float const ztop(water_plane_z - (dist_water ? 0.1 : 0.0));
			float const x1(bcube.d[0][dim ? 0 : dir]), x2(bcube.d[0][dim ? 1 : dir]), y1(bcube.d[1][dim ? dir : 0]), y2(bcube.d[1][dim ? dir : 1]);
			wverts.emplace_back(point(x1, y1, mzmin));
			wverts.emplace_back(point(x2, y2, mzmin));
			wverts.emplace_back(point(x2, y2, ztop));
			wverts.emplace_back(point(x1, y1, ztop));
		} // for dir
	} // for dim
	if (!wverts.empty()) {
		if (!textures_already_set && weight_tid > 0) {bind_textures();} // in the rare case where textures haven't been created, I guess we don't bind them
		int const loc(s.get_uniform_loc("htex_scale"));
		if (loc >= 0) {s.set_uniform_float(loc, 0.0);} // disable height texture
		draw_quad_verts_as_tris(wverts);
		if (loc >= 0) {s.set_uniform_float(loc, 1.0);} // enable height texture
	}
}


bool tile_t::is_water_visible() const {
	return (!is_distant && has_water() && rel_dist_to_camera_xy_lt(DRAW_DIST_TILES) && is_visible());
}
void tile_t::draw_water(shader_t &s, float z) const {

	if (!is_water_visible()) return;
	float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2)), xv2(xv1+(x2-x1)*deltax), yv2(yv1+(y2-y1)*deltay);
	bind_texture_tu(height_tid, 2);
	//bind_texture_tu(shadow_tid, 6); // Note: only needed if ENABLE_WATER_SHADOWS is enabled in the water plane shader
	bind_and_setup_shadow_map(s); // okay if shadow maps haven't been created yet
	draw_one_tquad(xv1, yv1, xv2, yv2, z, (use_water_plane_tess() ? GL_PATCHES : GL_TRIANGLE_FAN));
}


bool tile_t::check_sphere_collision(point &pos, float sradius, bool inc_dtrees, bool inc_ptrees, bool inc_scenery) const { // pos is in camera space

	if (is_distant || !contains_point(pos)) return 0;
	if (pos.z > get_tile_zmax() + sradius)  return 0; // sphere is completely above the tile
	bool coll(0);

	if (inc_ptrees && !pine_trees.empty()) {
		pos  -= ptree_off.get_xlate(); // Note: pos adj is required because pos is modified
		coll |= pine_trees.check_sphere_coll(pos, sradius);
		pos  += ptree_off.get_xlate();
	}
	if (inc_dtrees && !decid_trees.empty()) {
		pos  -= dtree_off.get_xlate();
		coll |= decid_trees.check_sphere_coll(pos, sradius);
		pos  += dtree_off.get_xlate();
	}
	if (inc_scenery && scenery.generated) {
		pos  -= scenery_off.get_xlate();
		coll |= scenery.check_sphere_coll(pos, sradius);
		pos  += scenery_off.get_xlate();
	}
	return coll;
}

bool tile_t::check_cube_int_trees(cube_t const &c) const { // cube is in camera space; deciduous trees only for now
	if (decid_trees.empty()) return 0;
	return decid_trees.check_cube_int(c - dtree_off.get_xlate());
}


int tile_t::get_tid_under_point(point const &pos) const {

	if (is_distant || !contains_point(pos)) return -1;
	int const xpos(max(0, min((int)size, (get_xpos(pos.x) - x1 - xoff + xoff2)))); // min/max not needed?
	int const ypos(max(0, min((int)size, (get_ypos(pos.y) - y1 - yoff + yoff2))));
	unsigned const ix(4*(ypos*stride + xpos));
	assert(ix < weight_data.size());
	unsigned max_weight(0), max_weight_ix(0), weight_sum(0);

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const w(weight_data[ix + i]); // mesh weight + tree dirt
		weight_sum += w;
		if (w > max_weight) {max_weight = w; max_weight_ix = i;}
	}
	if ((255 - weight_sum) > max_weight) {max_weight_ix = 4;} // 5th weight (255 - sum(w)) is max
	return lttex_dirt[max_weight_ix].id;
}


bool tile_t::line_intersect_mesh(point const &v1, point const &v2, float &t, int &xpos, int &ypos, float inc_trees) const {

	if (is_distant) return 0; // Note: this can be made to work, but won't work as-is
	assert(!is_nan(v1) && !is_nan(v2));
	point v1c(v1), v2c(v2); // clipped verts
	
	if (do_line_clip(v1c, v2c, get_mesh_bcube().d)) {
		// similar to mesh_intersector::line_intersect_surface_fast()
		int const xp1(get_xpos(v1c.x) - x1 - xoff + xoff2), yp1(get_ypos(v1c.y) - y1 - yoff + yoff2);
		int const xp2(get_xpos(v2c.x) - x1 - xoff + xoff2), yp2(get_ypos(v2c.y) - y1 - yoff + yoff2);
		int const dx(xp2 - xp1), dy(yp2 - yp1), steps(max(1, max(abs(dx), abs(dy))));
		assert(steps < 10000); // sanity check
		double const dz((double)v2c.z - (double)v1c.z), xinc(dx/(double)steps), yinc(dy/(double)steps), zinc(dz/(double)steps);
		double x(xp1), y(yp1), z(v1c.z - 0.1*fabs(zinc)); // z offset required to avoid problems with zval at bcube.z1

		for (int k = 0; k <= steps; ++k) {
			int const ix((int)x), iy((int)y);

			if (ix >= 0 && iy >= 0 && ix <= (int)size && iy <= (int)size && zvals[iy*zvsize + ix] > z) {
				// Note: we use z instead of zvals here because zvals may be much too high if we enter this tile while the line is under the mesh
				float const cur_t(((z - 0.5*zinc) - v1.z)/(double(v2.z) - double(v1.z))); // t relative to original v1, v2

				if (cur_t >= 0.0 && cur_t <= 1.0) {
					xpos = x1 + ix;
					ypos = y1 + iy;
					t    = cur_t;
					return 1;
				}
			}
			x += xinc; y += yinc; z += zinc;
		}
	}
	if (inc_trees) {
		if (!pine_trees.empty()) {
			vector3d const xlate(ptree_off.get_xlate());
			if (pine_trees.line_intersect((v1 - xlate), (v2 - xlate), &t)) {return 1;}
		}
		//if (!decid_trees.empty()) {} // check decid trees with -= dtree_off.get_xlate()
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
		float const mzval(get_exact_zval(pos.x, pos.y)), zstop(max(mzval, water_plane_z));
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
		if (time > int(fticks*LITNING_TIME2)) {clear();} // check for end time
		else if (animate2) {time += iticks;}
	}
	else if (animate2 && (rnum <= 1 || (rand()%rnum) == 0)) {gen();}
}

void lightning_strike_t::draw() const {

	if (!enabled()) return;
	path.draw_lines(1); // uses a custom shader with no fog; fade_ends=1
	colorRGBA const ambient(path.color*0.2);
	float const radius(0.4*get_scaled_tile_radius());
	set_colors_and_enable_light(LIGHTNING_LIGHT, ambient, path.color);
	setup_gl_light_atten(LIGHTNING_LIGHT, 0.1, 0.0, 1.0f/(radius*radius));
	set_gl_light_pos(LIGHTNING_LIGHT, get_pos(), 1.0); // point light source position
	tt_lightning_enabled = 1;
}

void lightning_strike_t::end_draw() const {

	disable_light(LIGHTNING_LIGHT); // even if not currently enabled, in case it was enabled before an update
	tt_lightning_enabled = 0;
}


// *** tile_draw_t ***


tile_draw_t::tile_draw_t() : lod_renderer(USE_TREE_BILLBOARDS) {
	assert(MESH_X_SIZE == MESH_Y_SIZE && X_SCENE_SIZE == Y_SCENE_SIZE);
}

void tile_draw_t::clear(bool no_regen_buildings) {

	clear_vbos_tids(); // needed to clear vbo, ivbo, and free list
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear();} // may not be necessary
	to_draw.clear();
	tiles.clear();
	shadow_recomp_queue.clear();
	if (!no_regen_buildings && !have_cities()) {buildings_valid = 0;} // can't regenerate buildings after cities and cars have been placed
}

void tile_draw_t::insert_tile(tile_t *tile) {
	bool const did_ins(tiles.insert(make_pair(tile->get_tile_xy_pair(), tile)).second);
	assert(did_ins);
}

void tile_draw_t::free_compute_shader() {
	for (auto i = height_gens.begin(); i != height_gens.end(); ++i) {i->clear_context();}
}

float tile_draw_t::update(float &min_camera_dist) { // view-independent updates; returns terrain zmin

	//highres_timer_t timer("TT Update");
	unsigned const max_tile_gen_per_frame = 16; // higher = less overall gen time (more parallel), but longer wait for first render
	unsigned const max_cpu_tiles          = 3; // 0 = GPU only
	unsigned const max_defer_tiles        = 8; // 0 = disable
	if (height_gens.empty()) {height_gens.resize(max(max_defer_tiles, 1U));}

	if (terrain_hmap_manager.maybe_load(mh_filename_tt, (invert_mh_image != 0))) {
		read_default_hmap_modmap();
		force_onto_surface_mesh(surface_pos); // move camera onto newly loaded terrain so that the first drawn frame is correct
	}
	else if (tiled_terrain_gen_heightmap_sz > 0) {
		terrain_hmap_manager.proc_gen_heightmap(tiled_terrain_gen_heightmap_sz);
		read_default_hmap_modmap();
		// since the heightmap values should be the same as single point queries, we don't need to re-calculate the player's zval
	}
	if (!buildings_valid) {
		gen_buildings();
		gen_city_details(); // after building generation
		buildings_valid = 1;
	}
	auto_calc_model_zvals(); // must be done after heightmap loading but before any tiles are created
	to_draw.clear();
	terrain_zmin = FAR_DISTANCE;
	grass_tile_manager.update(); // every frame, even if not in tiled terrain mode?
	assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
	point const cpos(get_camera_pos()), camera(cpos - get_tiled_terrain_model_xlate());
	int const tile_radius(int(CREATE_DIST_TILES*TILE_RADIUS) + 1);
	int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
	int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
	int const x2( tile_radius + toffx), y2( tile_radius + toffy);
	unsigned const init_tiles((unsigned)tiles.size());
	bool const create_buildings_first(FLATTEN_BUILDING_TILE && using_tiled_terrain_hmap_tex());
	unsigned num_erased(0);
	min_camera_dist = FAR_DISTANCE;
	// Note: we may want to calculate distant low-res or larger tiles when the camera is high above the mesh

	if (!to_gen_zvals.empty()) {
		assert(to_gen_zvals.size() <= height_gens.size());

		for (unsigned i = 0; i < to_gen_zvals.size(); ++i) { // tiles were waiting on zval generation (async)
			tile_t *tile(to_gen_zvals[i].second);
			tile->create_zvals(height_gens[i], 0); // wait for zvals to be generated
			insert_tile(tile); // zvals have been generated
		}
		to_gen_zvals.clear();
	}
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ) { // update tiles and free old tiles (Note: no ++i)
		if (!i->second->update_range(smap_manager)) { // delete this tile
			remove_buildings_tile(i->first.x, i->first.y); // required to avoid memory leak when player teleports to a new location
			i->second->clear();
			tiles.erase(i++);
			++num_erased;
		} else {++i;}
	}
	for (int y = y1; y <= y2; ++y ) { // create new tiles
		for (int x = x1; x <= x2; ++x ) {
			tile_xy_pair const txy(x, y);
			if (tiles.find(txy) != tiles.end()) continue; // already exists
			tile_t tile(get_tile_size(), x, y);
			if (!tile.rel_dist_to_camera_xy_lt(CREATE_DIST_TILES)) continue; // too far away to create
			tile_t *new_tile(new tile_t(tile));
			to_gen_zvals.push_back(make_pair(new_tile->get_draw_priority(), new_tile));
			// in this mode, we need to place buildings and flatten the heightmap before calculating tile heights
			if (create_buildings_first) {create_buildings_tile(x, y, 1);}
		} // for x
	} // for y
	//if (to_gen_zvals.size() < max_cpu_tiles) {to_gen_zvals.clear();} // block until at least max_cpu_tiles tiles to generate (lower average gen time, but causes more slow frames/lag)
	unsigned const num_to_gen(to_gen_zvals.size());
	unsigned gen_this_frame(min(num_to_gen, max_tile_gen_per_frame));
	bool const gpu_mode(mesh_gen_mode >= MGEN_SIMPLEX_GPU);
	
	// to balance tile gen time across frames, generate a number of tiles equal to the average of this frame and the previous frame
	if (gen_this_frame > 1 && gen_this_frame < max_tile_gen_per_frame && inf_terrain_fire_mode == FM_NONE) { // disable this mode when editing mesh height to prevent visual artifacts
		gen_this_frame = min(gen_this_frame, (gen_this_frame + tiles_gen_prev_frame + 1)/2); // round up
	}
	tiles_gen_prev_frame = num_to_gen;

	if (num_to_gen == 0) {
		// do nothing
	}
	else if (gpu_mode && num_to_gen <= max_defer_tiles) { // async generation mode - delay until next frame
		unsigned tgz_pos(0);

		for (unsigned i = 0; i < num_to_gen; ++i) {
			tile_t *tile(to_gen_zvals[i].second);
			if (tile->create_zvals(height_gens[i], 1)) {insert_tile(tile);} // no_wait=1; zvals have been generated, insert tile and remove from to_gen_zvals
			else {to_gen_zvals[i] = to_gen_zvals[tgz_pos++];} // zvals are not ready, leave in to_gen_zvals and try again during the next update
		}
		to_gen_zvals.resize(tgz_pos);
	}
	else {
		// if there are fewer than 4 tiles to generate, use CPU simplex rather than GPU simplex to avoid stalling/flusing the graphics pipeline
		int const prev_mesh_gen_mode(mesh_gen_mode);
		if (gpu_mode && gen_this_frame <= max_cpu_tiles) {mesh_gen_mode = MGEN_SIMPLEX;} // GPU simplex => CPU simplex
		if (gen_this_frame < num_to_gen) {sort(to_gen_zvals.begin(), to_gen_zvals.end());} // sort by priority if not all generated

		for (unsigned i = 0; i < num_to_gen; ++i) {
			tile_t *tile(to_gen_zvals[i].second);
			if (i >= gen_this_frame) {delete tile; continue;} // delete these tiles - they will be created in a later frame
			tile->create_zvals(height_gens[0], 0); // generate these tiles
			insert_tile(tile);
		}
		to_gen_zvals.clear();
		mesh_gen_mode = prev_mesh_gen_mode;
	}
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) { // calculate terrain_zmin and updated building tiles
		float const rel_dist(i->second->get_rel_dist_to_camera());

		if (rel_dist <= DRAW_DIST_TILES) {
			min_eq(terrain_zmin, i->second->get_zmin());
			if (!camera_surf_collide) {min_camera_dist = min(min_camera_dist, i->second->get_min_dist_to_pt(cpos, 0, 0));}
			create_buildings_tile(i->first.x, i->first.y, 0); // create, or re-create if create_buildings_first; should already be flat
		}
		else if (rel_dist > CLEAR_DIST_TILES) {remove_buildings_tile(i->first.x, i->first.y);}
	}
	if (DEBUG_TILES && (tiles.size() != init_tiles || num_erased > 0)) {
		cout << "update: tiles: " << init_tiles << " to " << tiles.size() << ", erased: " << num_erased << endl;
	}
	// Note: could skip shadow computation (but not weight calc/texture upload) if (max(sun_pos.z,  last_sun.z) > zbottom) or (sun.get_norm().z > 0.9) or something like that
	static point last_sun(all_zeros), last_moon(all_zeros);
	bool sun_change (sun_pos  != last_sun  && light_factor >= 0.4);
	bool moon_change(moon_pos != last_moon && light_factor <= 0.6 && max(moon_pos.z, last_moon.z) > zmin); // only when the moon is up
	float const toler = 0.9999; // only update when sun/moon has changed significantly
	sun_change  &= (dot_product(sun_pos.get_norm(),  last_sun.get_norm())  < toler);
	moon_change &= (dot_product(moon_pos.get_norm(), last_moon.get_norm()) < toler);

	if (mesh_shadows_enabled() && (sun_change || moon_change) && shadow_recomp_queue.empty()) { // light source change
		if (auto_time_adv && !moon_change) { // auto time advance shadow map update for sun change only - triger a shadow recompute
			for (auto i = tiles.begin(); i != tiles.end(); ++i) { // triger a shadow recompute
				shadow_recomp_queue.emplace_back(-p2p_dist(sun_pos, i->second->get_center()), i->second->get_tile_xy_pair());
			}
			sort(shadow_recomp_queue.begin(), shadow_recomp_queue.end()); // sort by decreasing distance to light source
		}
		else { // invalidate and recompute all shadows on moon change (infrequent) or user sun pos change
			for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear_shadows(sun_change, moon_change);}
		}
		last_sun  = sun_pos;
		last_moon = moon_pos;
	}
	unsigned num_shadow_updates = 12; // hard max per frame

	while (!shadow_recomp_queue.empty() && num_shadow_updates > 0) { // perform some queued shadow map updates, starting at light source
		tile_xy_pair const tp(shadow_recomp_queue.back().second);
		shadow_recomp_queue.pop_back();
		tile_map::const_iterator it(tiles.find(tp));
		if (it == tiles.end()) continue; // tile no longer exists/was deleted
		// recompute shadows; tiles feeding in (closer to the light) should have already been calculated
		it->second->clear_shadows(1, 0); // update sun shadows only
		it->second->check_shadow_map_and_normal_texture(1); // no_push=1
		--num_shadow_updates;
	}
	// Note: we could regen trees and scenery if water was just turned on to remove underwater vegetation
	//if ((GET_TIME_MS() - timer1) > 100) {PRINT_TIME("Tiled Terrain Update");}
	return terrain_zmin;
}

float tile_draw_t::get_actual_zmin() const {return min(zmin, terrain_zmin);}


float const mesh_tex_cscale [NTEX_DIRT] = {1.0, 1.0, TT_GRASS_COLOR_SCALE, 0.5, 1.0}; // darker grass and rock
float const mesh_tex_scale  [NTEX_DIRT] = {1.0, 1.0, 4.0, 1.0,  1.0};
int const normal_tids_dirt  [NTEX_DIRT] = {ROCK2_NORMAL_TEX, ROCK3_NORMAL_TEX, /*DIRT_NORMAL_TEX*/ROCK3_NORMAL_TEX, ROCK1_NORMAL_TEX, ROCK_NORMAL_TEX};
bool const disable_for_grass[NTEX_DIRT] = {0, 0, 1, 0, 0};

colorRGBA get_avg_color_for_landscape_tex(unsigned id) {
	assert(id < NTEX_DIRT);
	return texture_color(lttex_dirt[id].id) * mesh_tex_cscale[id];
}

/*static*/ void tile_draw_t::setup_terrain_textures(shader_t &s, unsigned start_tu_id) {

	unsigned const base_tsize(NORM_TEXELS);

	for (int i = 0; i < NTEX_DIRT; ++i) {
		int const tid(lttex_dirt[i].id);
		int const nm_tid((disable_for_grass[i] && is_grass_enabled()) ? FLAT_NMAP_TEX : normal_tids_dirt[i]);
		float const tscale(mesh_tex_scale[i]*float(base_tsize)/float(get_texture_size(tid, 0))); // assumes textures are square
		unsigned const tu_id(start_tu_id + i), nm_tu_id(i + 16); // tu_id 16-20 for normal maps
		select_texture(tid, tu_id);
		select_texture(nm_tid, nm_tu_id);
		assert(tu_id <= 9); // must map to a single character
		string tu_id_str;
		tu_id_str.push_back('0' + tu_id);
		s.add_uniform_int(  ("tex"  + tu_id_str).c_str(), tu_id);
		s.add_uniform_float(("ts"   + tu_id_str).c_str(), tscale);
		s.add_uniform_float(("cs"   + tu_id_str).c_str(), mesh_tex_cscale[i]);
		s.add_uniform_int(("nm_tex" + tu_id_str).c_str(), nm_tu_id);
	} // for i
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
	if (cloud_shadows) {s.set_prefix("#define APPLY_CLOUD_SHADOWS", lighting_shader);}
}

float get_cloud_coverage(float cloud_cover_factor) {return (is_cloudy ? 1.0 : (cloud_cover + cloud_cover_factor*(1.0 - cloud_cover)));}

void setup_cloud_plane_uniforms(shader_t &s, float cloud_cover_factor=0.535, bool match_cloud_layer=0) {

	float cloud_zmax;
	if (match_cloud_layer) {cloud_zmax = get_cloud_zmax();} // follows the camera zval - matches the drawn cloud layer but moves clouds on the terrain
	else {cloud_zmax = 0.5f*(zmin + zmax) + max(zmax, CLOUD_CEILING);} // fixed z value - independent of camera z so stays in place, but disagrees with drawn clouds
	set_cloud_uniforms(s, 9);
	s.add_uniform_float("cloud_scale",   get_cloud_coverage(cloud_cover_factor));
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

void setup_tile_shader_shadow_map(shader_t &s) {

	for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) {
		string sm_tex_str("sm_tex");
		s.add_uniform_int(append_ix(sm_tex_str, i, 0), TILE_SMAP_START_TU_ID+i);
	}
	s.add_uniform_float("smap_atten_cutoff", get_smap_atten_val());
	s.add_uniform_float("z_bias", DEF_Z_BIAS);
}

void set_smap_enable_for_shader(shader_t &s, bool enable_smap, int shader_type) {
	s.set_prefix(make_shader_bool_prefix("use_shadow_map", enable_smap), shader_type);
	if (enable_smap) {s.set_prefix("#define ENABLE_SHADOW_MAP", shader_type);} // need to set this as well
}


// uses texture units 0-11 and 15 (12 if using hmap texture, 13-14 if using shadow maps) + 16-20 for normal maps
// Note: could be static, except uses get_actual_zmin()
void tile_draw_t::setup_mesh_draw_shaders(shader_t &s, bool reflection_pass, bool enable_shadow_map) const {

	bool const has_water(is_water_enabled() && !reflection_pass);
	lighting_with_cloud_shadows_setup(s, 1, (cloud_shadows_enabled() && !reflection_pass));
	bool const water_caustics(has_water && !(display_mode & 0x80) && water_params.alpha < 1.5);
	bool const use_normal_map(!reflection_pass); // enabled by default
	bool const rain_mode((precip_mode > 0) && temperature > W_FREEZE_POINT);
	bool const triplanar_tex = tt_triplanar_tex; // slower, looks somewhat better on steep terrain
	if (has_water      ) {s.set_prefix("#define HAS_WATER",            1);} // FS
	if (water_caustics ) {s.set_prefix("#define WATER_CAUSTICS",       1);} // FS
	if (use_normal_map ) {s.set_prefix("#define USE_NORMAL_MAP",       1);} // FS
	if (reflection_pass) {s.set_prefix("#define REFLECTION_MODE",      1);} // FS
	if (triplanar_tex  ) {s.set_prefix("#define TRIPLANAR_TEXTURE",    1);} // FS
	if (combined_gu)     {s.set_prefix("#define ADD_INDIR_REFL_LIGHT", 1);} // FS - add extra indir reflected light in combined_gu mode to offset the low ambient
	set_smap_enable_for_shader(s, enable_shadow_map, 1); // FS
	s.set_vert_shader("tiled_mesh");
	s.set_frag_shader("water_fog.part*+linear_fog.part+perlin_clouds.part*+ads_lighting.part*+shadow_map.part*+detail_normal_map.part+triplanar_texture.part+tiled_mesh");
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("weights_tex", 0);
	s.add_uniform_int("detail_tex",  1);
	s.add_uniform_int("normal_tex",  7);
	s.add_uniform_int("shadow_tex",  15);
	s.add_uniform_float("normal_z_scale", (reflection_pass ? -1.0 : 1.0));
	s.add_uniform_float("spec_offset", (rain_mode ? 0.5 : 0.0)); // increase specular during rain
	set_noise_tex(s, 8);
	setup_cloud_plane_uniforms(s);
	setup_terrain_textures(s, 2); // TU_ID=2-7
	set_tile_xy_vals(s);
	s.add_uniform_int("height_tex", 12);
	if (use_normal_map   ) {setup_detail_normal_map(s, 2.0);}
	if (enable_shadow_map) {setup_tile_shader_shadow_map(s);}

	if (has_water) {
		set_water_plane_uniforms(s);
		setup_shader_underwater_atten(s, WATER_COL_ATTEN*mesh_scale);

		if (water_caustics) {
			select_texture(WATER_CAUSTIC_TEX, 10);
			s.add_uniform_int("caustic_tex", 10);
			s.add_uniform_float("caustics_weight", (1.5 - water_params.alpha));
		}
		setup_water_plane_texgen(8.0, 2.5, s, 2); // increase texture scale and change AR since the caustics texture is sparser than the water texture
	}
	else { // or just disable water fog calculation in the vertex shader (water_fog.part)?
		s.add_uniform_float("water_plane_z", (reflection_pass ? water_plane_z : get_actual_zmin())); // used for fog calculation/clipping
	}
	set_landscape_texgen(1.0, -MESH_X_SIZE/2, -MESH_Y_SIZE/2, MESH_X_SIZE, MESH_Y_SIZE, s, 1);
	s.add_uniform_float("triplanar_texture_scale", 1.0f/(X_SCENE_SIZE + Y_SCENE_SIZE));
}


bool tile_draw_t::can_have_reflection_recur(tile_t const *const tile, point const corners[3], unsigned dim_ix) {

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
	if (tile->vis_ref_call) return 0; // already seen
	tile->vis_ref_call = 1;
	bool ret(0);
		
	for (unsigned d = 0; d < 2; ++d) {
		int delta[2] = {0,0};
		if      (camera[d] < bcube.d[d][0]) {delta[d] = -1;}
		else if (camera[d] > bcube.d[d][1]) {delta[d] =  1;}
		else                                {continue;}
		tile_t const *const adj(get_tile_from_xy(tile->get_tile_xy_pair(delta[0], delta[1])));
		if (!adj || !adj->is_visible()) continue;
		if (can_have_reflection_recur(adj, corners, d)) {ret = 1; break;}
	}
	tile->vis_ref_call = 0;
	return ret;
}


bool tile_draw_t::can_have_reflection(tile_t const *const tile) {

	if (!is_water_enabled()) return 0;
	if (tile->all_water  ()) return 0;
	if (tile->has_water  ()) return 1;
	// return 1 if the camera's tile contains water?
	cube_t bcube(tile->get_bcube());
	point const camera(get_camera_pos()), center(bcube.get_cube_center());
	point corners[3]; // {x, y, closest}

	for (unsigned d = 0; d < 2; ++d) {
		corners[ d][d] = bcube.d[d][center[d] < camera[d]]; // 0,0
		corners[!d][d] = bcube.d[d][center[d] > camera[d]]; // 1,1
		corners[ d][2] = bcube.d[2][0];
	}
	corners[2]   = bcube.closest_pt(camera);
	corners[2].z = tile->get_zmax();
	return can_have_reflection_recur(tile, corners, 2);
}


unsigned get_building_models_gpu_mem();
unsigned get_dlights_smap_gpu_mem();
unsigned get_vbo_ring_buffers_size();
unsigned get_quad_ix_buffer_size();

uint64_t tile_draw_t::show_debug_stats(bool calc_mem_only) const {
	unsigned num_trees(0), num_smaps(0);
	uint64_t mem(0), tree_mem(0), smap_mem(0);

	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		tile_t *const tile(i->second.get());
		mem       += tile->get_gpu_mem (); // Note: includes smap_mem
		tree_mem  += tile->get_tree_mem();
		smap_mem  += tile->get_smap_mem();
		num_smaps += tile->count_shadow_maps();
		num_trees += tile->num_pine_trees() + tile->num_decid_trees();
	}
	unsigned const dtree_mem(tree_data_manager.get_gpu_mem()), ptree_mem(get_tree_inst_gpu_mem()), grass_mem(grass_tile_manager.get_gpu_mem());
	unsigned const smap_free_list_mem(smap_manager.get_free_list_mem_usage()), dlights_smap_mem(get_dlights_smap_gpu_mem());
	unsigned const texture_mem( get_loaded_textures_gpu_mem());
	unsigned const building_mem(get_buildings_gpu_mem_usage());
	unsigned const models_mem(get_city_model_gpu_mem() + get_loaded_models_gpu_mem() + get_building_models_gpu_mem());
	unsigned const frame_buf_mem(13*window_width*window_height); // RGB8 (as 32 bits?) front buffer + RGB8 back buffer + 32-bit depth buffer + 8 bit stencil buffer
	unsigned const ring_buf_mem(get_vbo_ring_buffers_size()), quad_ix_buf_mem(get_quad_ix_buffer_size()); // these should be small (8/16MB, ~5MB)

	if (vbo) {
		unsigned const tile_size(get_tile_size());
		mem += 2ULL*tile_size*(tile_size+1ULL)*sizeof(point);
		for (unsigned i = 0; i < NUM_LODS; ++i) {mem += 4ULL*(tile_size>>i)*(tile_size>>i)*sizeof(unsigned);} // approximate
	}
	uint64_t const tot_mem(mem + dtree_mem + ptree_mem + grass_mem + smap_free_list_mem + dlights_smap_mem + texture_mem + building_mem + models_mem +
		frame_buf_mem + room_geom_mem + ring_buf_mem + quad_ix_buf_mem);
	if (calc_mem_only) return tot_mem;

	cout << "tiles drawn: " << to_draw.size() << " of " << tiles.size() << ", trees drawn: " << num_trees << ", shadow maps: " << num_smaps
		<< ", GPU MB: " << in_mb(tot_mem)
		<< ", tile MB: " << in_mb(mem - smap_mem) << ", tree CPU MB: " << in_mb(tree_mem) << ", tree GPU MB: " << in_mb((unsigned long long)dtree_mem + ptree_mem)
		<< ", grass MB: " << in_mb(grass_mem) << ", smap MB: " << in_mb(smap_mem) << ", smap free list MB: " << in_mb(smap_free_list_mem)
		<< ", dlights smap mem MB: " << in_mb(dlights_smap_mem) << ", frame buf MB: " << in_mb(frame_buf_mem) << ", texture MB: " << in_mb(texture_mem)
		<< ", building MB: " << in_mb(building_mem) << ", room_geom MB: " << in_mb(room_geom_mem) << ", model MB: " << in_mb(models_mem) << endl;
	//show_gpu_mem_info(); // shows total and available video memory
	return tot_mem;
}


void tile_draw_t::pre_draw() { // view-dependent updates/GPU uploads

	//timer_t timer("TT Pre-Draw");
	vector<tile_t *> to_update, to_update_shadows, to_gen_trees;
	assert((vbo == 0) == (ivbo == 0)); // either neither or both are valid
	get_empty_smap_tid(); // we're going to need this later, so make sure to allocate the texture first so that it doesn't invalidate TU 0 mid-tile draw

	// handle clearing of tile shadow maps
	vector<clear_area_t> to_clear_next_frame;

	for (clear_area_t const &i : tile_smaps_to_clear) {
		invalidate_tile_smap_in_region(i);
		if (i.clear_next_frame) {to_clear_next_frame.emplace_back(i, 0);} // insert again with clear_next_frame=0
	}
	tile_smaps_to_clear = to_clear_next_frame;
	
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
		if (!tile->rel_dist_to_camera_xy_lt(DRAW_DIST_TILES)) continue; // too far to draw
		//if (display_mode & 0x20) {tile->clear_shadow_map(&smap_manager);} // useful for perf testing
		bool const is_visible(tile->is_visible());

		if (shadow_map_enabled()) {
			if (is_visible || tile->is_smap_bounds_visible()) {to_update_shadows.push_back(tile);} // Note: using current camera view frustum
			else {tile->setup_shadow_maps(smap_manager, 1);} // cleanup_only=1 (only clear shadow maps to increase LOD levels)
		}
		if (is_visible) {
			if (tile->can_have_trees()) { // no trees in water or distant tiles
				if (tile->can_have_pine_palm_trees() && !tile->pine_trees_generated()) {to_gen_trees.push_back(tile);}
				if (decid_trees_enabled()) {tile->gen_decid_trees_if_needed();}
			}
			to_update.push_back(tile);
		}
	} // for i
	if (enable_instanced_pine_trees() && !to_gen_trees.empty()) {create_pine_tree_instances();}
	//RESET_TIME;
	// don't use parallel tree gen for a single tile
#pragma omp parallel for schedule(dynamic,1) if (to_gen_trees.size() > 1)
	for (int i = 0; i < (int)to_gen_trees.size(); ++i) {to_gen_trees[i]->init_pine_tree_draw();}
	//if (!to_gen_trees.empty()) {PRINT_TIME("Gen Trees2");}
	assert(!height_gens.empty());
	
	for (vector<tile_t *>::iterator i = to_update.begin(); i != to_update.end(); ++i) {
		(*i)->pre_draw(height_gens[0]);

		if ((*i)->can_have_trees()) {
			(*i)->update_pine_tree_state(1);
			(*i)->update_decid_trees();
		}
		(*i)->update_scenery();
	} // for i
	for (vector<tile_t *>::iterator i = to_update_shadows.begin(); i != to_update_shadows.end(); ++i) { // after everything has been setup
		(*i)->setup_shadow_maps(smap_manager, 0); // cleanup_only=0
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

tile_draw_t::occluder_cubes_t::occluder_cubes_t(tile_t const *const tile_) : tile(tile_), bcube(tile->get_mesh_bcube()) {
	for (unsigned s = 0; s < 16; ++s) {
		cube_t &sub_cube(sub_cubes[s]);
		sub_cube = tile->get_mesh_sub_bcube((s>>2), (s&3));
		if (!camera_pdu.cube_visible(sub_cube)) {sub_cube.set_to_zeros(); continue;} // if cube is behind the player, it can't be an occluder; doesn't help much
		sub_cube.z2() = sub_cube.z1(); // cube below the bcube
		sub_cube.z1() = zmin;
	}
}

void tile_draw_t::draw(int reflection_pass) { // reflection_pass: 0=none, 1=water plane Z, 2=building mirror

	if (player_cant_see_outside_building()) return; // no need to draw tiles if player in extended basement or parking garage
	//timer_t timer("TT Draw");
	shadow_maps_disabled = 0; // reset for this frame
	to_draw.clear();
	occluded_tiles.clear();

	// determine potential occluders
	point const camera(get_camera_pos());
	occluder_cubes.clear();
	if (!building_occluder.is_all_zeros()) {occluder_cubes.push_back(building_occluder);}

	if ((display_mode & 0x08) && (display_mode & 0x01) && check_tt_mesh_occlusion) { // check occlusion when occlusion culling and mesh are enabled
		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
			tile_t *const tile(i->second.get());
			if (tile->use_as_occluder()) {occluders.emplace_back(tile);}
		}
	}
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		tile_t *const tile(i->second.get());
		float const dist(tile->get_rel_dist_to_camera());
		if (dist > DRAW_DIST_TILES || !tile->is_visible()) continue;
		if (tile->was_last_occluded()) continue; // occluded in the shadow pass
		if (reflection_pass == 1 && !can_have_reflection(tile)) continue; // check for water plane Z reflections only
		if (!occluder_cubes.empty() && check_cube_occluded(tile->get_bcube(), occluder_cubes, camera)) continue; // building wall/floor/ceiling occlusion

		if (!occluders.empty() && !tile->was_last_unoccluded()) {
			occluder_pts_t tile_os, sub_tile_os;
			tile_os.calc_cube_top_points(tile->get_bcube());
			bool tile_occluded(1);
			occluder_ixs.clear();

			for (auto j = occluders.begin(); j != occluders.end(); ++j) {
				if (j->tile == tile) continue; // no self-occlusion
				bool intersected(0); // will always be true for the camera's current tile
				for (unsigned d = 0; d < 4 && !intersected; ++d) {intersected |= check_line_clip(tile_os.cube_pts[d], camera, j->bcube.d);}
				if (intersected) {occluder_ixs.push_back(j - occluders.begin());}
			} // for j
			for (unsigned t = 0; t < 16; ++t) {
				cube_t const sub_cube(tile->get_sub_bcube((t>>2), (t&3)));
				sub_tile_os.calc_cube_top_points(sub_cube);
				point const center(sub_cube.get_cube_center());
				bool sub_tile_occluded(0);

				for (auto S = occluder_ixs.begin(); S != occluder_ixs.end() && !sub_tile_occluded; ++S) {
					occluder_cubes_t const &oc(occluders[*S]);
					if (!check_line_clip(center, camera, oc.bcube.d)) continue; // if not occluded by the full tile, then won't be occluded by a sub-cube of it

					for (unsigned s = 0; s < 16 && !sub_tile_occluded; ++s) {
						cube_t const &occluder(oc.sub_cubes[s]);
						if (occluder.is_all_zeros()) continue; // invalid cube, skip
						sub_tile_occluded = 1;

						for (unsigned d = 0; d < 4 && sub_tile_occluded; ++d) {
							sub_tile_occluded &= check_line_clip(sub_tile_os.cube_pts[d], camera, occluder.d);
						}
					} // for s
				} // for S
				if (!sub_tile_occluded) {tile_occluded = 0; break;} // Note: must mark as not occluded, even if this sub-tile is not visible to the camera
			} // for t
			tile->set_last_occluded(tile_occluded);
			if (tile_occluded) {occluded_tiles.push_back(tile); continue;}
		} // check_occlusion
		to_draw.emplace_back(dist, tile);
	} // for i
	//cout << TXT(occluders.size()) << TXT(tiles.size()) << TXT(to_draw.size()) << TXT(occluded_tiles.size()) << endl; // 5, 341, 63, 34
	occluders.clear();
	occluder_ixs.clear();
	building_occluder.set_to_zeros(); // was used, may be reset during the next frame
	sort(to_draw.begin(), to_draw.end()); // sort front to back to improve draw time through depth culling

	if (display_mode & 0x01) { // draw visible tiles
		if (shadow_map_enabled()) {draw_tiles(reflection_pass, 1);} // shadow map pass
		draw_tiles(reflection_pass, 0); // non-shadow map pass
	}
	if (!player_cant_see_outside_building()) { // trees/scenerg/grass not visible when player is in the extended basement, parking garage, or attic
		if (pine_trees_enabled ()) {draw_pine_trees (reflection_pass);}
		if (decid_trees_enabled()) {draw_decid_trees(reflection_pass);}
	
		if (player_in_basement < 2 && player_in_attic < 2) { // vegetation/scenery/animals not visible when player is fully inside the basement or windowless attic
			if (scenery_enabled    ()) {draw_scenery    (reflection_pass);}
			if (is_grass_enabled   ()) {draw_grass      (reflection_pass);}
			if (ENABLE_ANIMALS)        {draw_animals    (reflection_pass);}
		}
	}
	if (DEBUG_TILES) {show_debug_stats(0);} // calc_mem_only=0
	//if ((GET_TIME_MS() - timer1) > 100) {PRINT_TIME("Draw Tiled Terrain");}
}

void tile_draw_t::end_lightning() const {lightning_strike.end_draw();} // in case it was enabled


void tile_draw_t::draw_tiles(int reflection_pass, bool enable_shadow_map) const {

	//timer_t timer("TT Draw Tiles");
	shader_t s;
	setup_mesh_draw_shaders(s, reflection_pass, enable_shadow_map);
	s.add_uniform_float("spec_scale", 1.0);
	s.clear_specular(); // in case we failed to clear it somewhere ahead
	s.enable_vnct_atribs(1, 0, 0, 0);
	enable_blend(); // for fog transparency
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(PRIMITIVE_RESTART_IX);
	int shader_locs[2] = {-1, -1};

	for (unsigned i = 0; i < to_draw.size(); ++i) {
		if (to_draw[i].second->using_shadow_maps() != enable_shadow_map) continue; // draw in another pass
		to_draw[i].second->draw(s, *this, ivbo_ixs, crack_ibuf, reflection_pass, shader_locs);
	}
	if (!enable_shadow_map && !reflection_pass && is_water_enabled()) { // draw all water caps (required, even for occluded tiles)
		for (auto i = occluded_tiles.begin(); i != occluded_tiles.end(); ++i) {(*i)->draw_water_cap(s, 0);}
	}
	bind_vbo(0, 1); // unbind index buffer
	glDisable(GL_PRIMITIVE_RESTART);
	disable_blend();

	if ((display_mode & 0x01) && !enable_shadow_map && !reflection_pass && draw_distant_water() && water_plane_z > terrain_zmin) {
		bind_2d_texture(BLACK_TEX); // all snow? at least it's set to something valid
		bind_texture_tu(WHITE_TEX, 15); // shadow_map texture, use something determinsitic (not that it matters visually)
		select_texture(FLAT_NMAP_TEX, 7); // normal_map texture
		disable_shadow_maps(s);
		int const loc(s.get_uniform_loc("htex_scale"));
		if (loc >= 0) {s.set_uniform_float(loc, 0.0);} // disable height texture
		draw_distant_mesh_bottom(terrain_zmin); // Note: textures from last drawn tile are bound, but don't really affect the results
		if (loc >= 0) {s.set_uniform_float(loc, 1.0);} // enable height texture
	}
	s.end_shader();
	//if (display_mode & 0x10) {draw_smap_debug_vis();} // TESTING
}

vector4d vec4_from_cube_xy(cube_t const &c) {return vector4d(c.x1(), c.y1(), c.x2(), c.y2());}

void tile_draw_t::draw_tiles_shadow_pass(point const &lpos, tile_t const *const tile) { // not const because creates height_tid

	//timer_t timer("Draw Shadow Pass");
	assert(tile != nullptr); // must have a valid dest tile
	shader_t s;
	s.set_vert_shader("tiled_mesh_shadow");
	s.set_frag_shader("color_only");
	s.begin_shader();
	set_tile_xy_vals(s);
	s.add_uniform_int("height_tex", 12);
	s.add_uniform_float("delta_z", -0.5*grass_length); // move mesh down by half the grass length so that the bottoms of the grass blades aren't shadowed

	if (pond_max_depth > 0.0 && tile->is_inside_city()) { // exclude city area so that ponds aren't shadowed
		cube_t const city_bcube(get_city_bcube_overlapping(tile->get_mesh_bcube_global())); // should be nonzero
		s.add_uniform_vector4d("exclude_box", vec4_from_cube_xy(city_bcube + get_camera_coord_space_xlate())); // global to camera space
		s.add_uniform_float("exclude_dz", -2.0*pond_max_depth);
	}
	else { // no exclude
		s.add_uniform_vector4d("exclude_box", vector4d());
		s.add_uniform_float("exclude_dz", 0.0);
	}
	s.enable_vnct_atribs(1, 0, 0, 0);
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(PRIMITIVE_RESTART_IX);
	float const recv_dist_sq(p2p_dist_xy_sq(lpos, tile->get_center()));
	cube_t const shadow_bcube(tile->get_shadow_bcube());
	bool const inc_adj_smap(get_buildings_max_extent() != zero_vector || get_road_max_len().x > 0.0);
	int const xlate_loc(s.get_uniform_loc("translate"));
	assert(xlate_loc >= 0);

	for (unsigned i = 0; i < to_draw.size(); ++i) {
		tile_t *const t(to_draw[i].second);

		if (p2p_dist_xy_sq(lpos, t->get_center()) <= recv_dist_sq || // check if this tile is closer to the light than the recv tile
		   (inc_adj_smap && t->get_bcube().intersects_xy(shadow_bcube))) // check for tile overlap when buildings/cities are involved
		{
			t->draw_shadow_pass(s, *this, ivbo_ixs, xlate_loc);
		}
	}
	bind_vbo(0, 1); // unbind index buffer
	glDisable(GL_PRIMITIVE_RESTART);
	s.end_shader();
}


void tile_draw_t::draw_shadow_pass(point const &lpos, tile_t *tile, bool decid_trees_only) {

	if (decid_trees_only && !decid_trees_enabled()) return;
	// frame time when forcing shadow map recreation every frame: 7ms base, 5ms for models, 5ms for decid trees, 4ms for everything else, 18ms total
	//highres_timer_t timer("Draw Shadow Pass"); // 0.99ms for cities with no cars
	float const orig_near_plane(camera_pdu.near_);
	bool const orig_fog_enabled(fog_enabled);
	camera_pdu.near_ = 0.0; // move the near clipping plane to zero to prevent clipping of tiles that are between the light and the target but not in the shadow frustum
	fog_enabled      = 0; // optimization?
	to_draw.clear();
	select_texture(WHITE_TEX); // make sure TU 0 is valid and not left as a shadow map texture

	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) { // 0.03ms
		if (decid_trees_only && i->second->num_decid_trees() == 0) continue; // no decid trees to draw, don't bother checking visibility
		if (!i->second->is_visible()) continue;
		if (!decid_trees_only) {i->second->update_pine_tree_state(1, 1);} // force high detail trees
		//i->second->update_decid_trees(); // not legal
		to_draw.push_back(make_pair(0.0, i->second.get())); // distance is unused so set to 0.0
	}
	if (!enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);} // enable depth clamping so that shadow casters aren't clipped by the shadow frustum

	if (!decid_trees_only) {
		draw_tiles_shadow_pass(lpos, tile); // 0.03ms

		if (pine_trees_enabled()) {
			draw_pine_trees(0, 1);
			for (auto i = to_draw.begin(); i != to_draw.end(); ++i) {i->second->update_pine_tree_state(1, 0);} // reset detail; 0.09ms
		}
		if (scenery_enabled()) {draw_scenery(0, 1);} // 0.12ms
		render_models(1, 0, 3, get_tiled_terrain_model_xlate()); // both transparent and opaque; VFC should work here for models; 0.45ms
	}
	if (decid_trees_enabled()) {draw_decid_trees(0, 1);} // 0.27ms / 5ms frame time
	if (!enable_depth_clamp) {glDisable(GL_DEPTH_CLAMP);}
	fog_enabled      = orig_fog_enabled;
	camera_pdu.near_ = orig_near_plane;
}


void tile_draw_t::draw_smap_debug_vis() const {
	cout << "smap sizes: ";
	shader_t s;
	s.begin_color_only_shader();
	ensure_outlined_polygons();
	for (auto const &i : to_draw) {i.second->draw_smap_debug_vis(s);}
	s.end_shader();
	set_fill_mode();
	cout << endl;
}


void tile_draw_t::draw_water(shader_t &s, float zval) const {
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->draw_water(s, zval);}
}


colorRGBA get_color_scale(float mag=1.0, float cloud_cover_factor=0.0) {
	float const lcscale(mag*(cloud_shadows_enabled() ? (1.0 - 0.5*get_cloud_coverage(cloud_cover_factor)) : 1.0)); // cloud scale is nominally 0.75 with cloud_cover_factor=0.5
	return colorRGBA(lcscale, lcscale, lcscale, 1.0);
}


/*static*/ void tile_draw_t::set_noise_tex(shader_t &s, unsigned tu_id) {
	select_texture(DITHER_NOISE_TEX, tu_id);
	s.add_uniform_int("noise_tex", tu_id);
}
/*static*/ void tile_draw_t::set_tree_dither_noise_tex(shader_t &s, unsigned tu_id) {
	set_noise_tex(s, tu_id);
	s.add_uniform_float("noise_tex_size", get_texture_size(DITHER_NOISE_TEX, 0));
}

/*static*/ void tile_draw_t::set_pine_tree_shader(shader_t &s, string const &vs, bool use_texgen) {
	shared_shader_lighting_setup(s, 0); // VS
	s.set_vert_shader((use_texgen ? "ads_lighting.part*+texture_gen.part+" : "ads_lighting.part*+") + vs);
	s.set_frag_shader("linear_fog.part+noise_dither.part+pine_tree");
	set_pine_tree_shader_post(s);
}
/*static*/ void tile_draw_t::set_pine_tree_shader_post(shader_t &s) {
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("branch_tex", 0);
	s.add_uniform_float("min_alpha", 0.5);
	s.add_uniform_color("color_scale", get_color_scale());
	set_tree_dither_noise_tex(s, 1); // TU=1
	check_gl_error(302);
}


void tile_draw_t::draw_pine_tree_bl(shader_t &s, bool branches, bool near_leaves, bool far_leaves, bool shadow_pass, bool reflection_pass, bool enable_smap, int xlate_loc) {

	for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
		to_draw[i].second->draw_pine_trees(s, to_draw_trunk_pts, branches, near_leaves, far_leaves, shadow_pass, reflection_pass, enable_smap, xlate_loc);
	}
}


void tile_draw_t::draw_pine_trees(bool reflection_pass, bool shadow_pass) { // and palm trees

	// far leaves
	if (!shadow_pass) {
		enable_blend(); // for fog transparency
		shader_t s;
		s.set_geom_shader("pine_tree_billboard"); // point => 1 quad
		set_pine_tree_shader(s, "pine_tree_billboard_gs", 0); // doesn't need texture_gen.part
		//set_pine_tree_shader(s, "pine_tree_billboard_auto_orient");
		disable_shadow_maps(s); // bind a valid shadow map
		s.add_uniform_float("radius_scale", calc_tree_size());
		s.add_uniform_float("ambient_scale", 1.5);
		s.set_specular(0.2, 8.0);
		draw_pine_tree_bl(s, 0, 0, 1, shadow_pass, reflection_pass, 0, s.get_uniform_loc("xlate"));
		s.clear_specular();
		s.end_shader();
		disable_blend();
	}
	// near leaves
	shader_t s;
	float const wind_mag(get_plant_leaf_wind_mag(shadow_pass));
	if (wind_mag > 0.0) {s.set_prefix("#define ENABLE_WIND", 0);} // VS
	if (shadow_pass   ) {s.set_prefix("#define NO_NOISE",    1);} // FS
	if (enable_instanced_pine_trees()) {s.set_prefix("#define ENABLE_INSTANCING", 0);} // VS
	bool const enable_smap(!shadow_pass && shadow_map_enabled()), enable_leaf_smap(!(display_mode & 0x0200) && enable_smap);

	if (enable_leaf_smap) { // per-pixel lighting/shadows
		shared_shader_lighting_setup(s, 1); // FS
		set_smap_enable_for_shader(s, 1, 1); // FS - smap always enabled
		s.set_prefix("#define NO_SHADOW_PCF",    1); // FS
		s.set_prefix("#define TEST_BACK_FACING", 1); // FS
		s.set_vert_shader("texture_gen.part+leaf_wind.part+pine_tree_ppl");
		s.set_frag_shader("ads_lighting.part*+linear_fog.part+noise_dither.part+shadow_map.part*+tiled_shadow_map.part*+pine_tree_ppl");
		set_pine_tree_shader_post(s);
		setup_tile_shader_shadow_map(s);
	}
	else {
		set_pine_tree_shader(s, "leaf_wind.part+pine_tree");
	}
	setup_leaf_wind(s, wind_mag, 0);
	int xlate_loc(-1);
	
	if (enable_instanced_pine_trees()) {
		s.add_uniform_float("vertex_scale", calc_tree_size()); // default is 0.8
		xlate_loc = s.get_attrib_loc("xlate");
		enable_instancing_for_shader_loc(xlate_loc);
	}
	if (!shadow_pass) {s.set_specular(0.2, 8.0);}
	draw_pine_tree_bl(s, 0, 1, 0, shadow_pass, reflection_pass, enable_leaf_smap, xlate_loc);
	assert(to_draw_trunk_pts.empty());
	if (!shadow_pass) {s.clear_specular();}
	if (xlate_loc >= 0) {disable_instancing_for_shader_loc(xlate_loc);}
	s.end_shader();

	// near trunks
	bool const enable_dlights(!shadow_pass && !reflection_pass && enable_tree_dlights());
	tree_branch_shader_setup(s, enable_smap, 0, shadow_pass, enable_dlights); // enable_opacity=0
	s.add_uniform_float("tex_scale_t", 5.0);
	draw_pine_tree_bl(s, 1, 0, 0, shadow_pass, reflection_pass, enable_smap, -1); // branches
	s.add_uniform_float("tex_scale_t", 1.0);
	s.end_shader();

	// distant/far trunks
	if (!shadow_pass && !to_draw_trunk_pts.empty()) { // color/texture already set above
		enable_blend(); // for fog transparency
		set_pine_tree_shader(s, "xy_billboard");
		select_texture(WHITE_TEX);
		s.set_cur_color(get_tree_trunk_color(T_PINE, 1));
		for (auto i = to_draw_trunk_pts.begin(); i != to_draw_trunk_pts.end(); ++i) {(*i)->draw_trunk_pts(s);}
		s.add_uniform_vector3d("xlate", zero_vector);
		s.end_shader();
		disable_blend();
	}
	to_draw_trunk_pts.clear();
}


void tile_draw_t::draw_decid_tree_bl(shader_t &s, tree_lod_render_t &lod_renderer, bool branches, bool leaves, bool reflection_pass, bool shadow_pass, bool enable_smap) {

	for (unsigned i = 0; i < to_draw.size(); ++i) { // near leaves
		to_draw[i].second->draw_decid_trees(s, lod_renderer, branches, leaves, reflection_pass, shadow_pass, enable_smap);
	}
}


void tile_draw_t::billboard_tree_shader_setup(shader_t &s) {

	shared_shader_lighting_setup(s, 1);
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("normal_map", 1);
	s.add_uniform_int("color_map",   0);
	s.add_uniform_int("tc_start_ix", 0);
	set_tree_dither_noise_tex(s, 2); // TU=2
}


void tile_draw_t::tree_branch_shader_setup(shader_t &s, bool enable_shadow_maps, bool enable_opacity, bool shadow_only, bool enable_dlights) {

	cube_t lights_bcube(all_zeros);

	if (enable_dlights) {
		lights_bcube = get_city_lights_bcube();
		if (lights_bcube.is_all_zeros()) {enable_dlights = 0;}
	}
	if (enable_opacity) {s.set_prefix("#define ENABLE_OPACITY",        1);} // FS
	if (enable_dlights) {s.set_prefix("#define ENABLE_DYNAMIC_LIGHTS", 1);} // FS
	if (rotate_trees && enable_dlights) {s.set_prefix("#define ENABLE_ROTATIONS", 0);} // VS
	s.setup_enabled_lights(3, 2); // FS; sun, moon, and lightning
	set_dlights_booleans(s, enable_dlights, 1, 1); // no_dl_smap=1
	if (!shadow_only) {setup_tt_fog_pre(s);}
	set_smap_enable_for_shader(s, enable_shadow_maps, 1); // FS
	s.set_vert_shader(enable_dlights ? "world_space_offset_rot.part+tiled_tree_branches" : "per_pixel_lighting");
	s.set_frag_shader(string("linear_fog.part+ads_lighting.part*+") + (enable_dlights ? "dynamic_lighting.part*+" : "") +
		"noise_dither.part+shadow_map.part*+tiled_shadow_map.part*+tiled_tree_branches");
	s.begin_shader();
	s.add_uniform_int("tex0", 0);
	//s.add_uniform_int("shadow_tex", 6);

	if (!shadow_only) {
		setup_tt_fog_post(s);
		s.add_uniform_color("ambient_tint", colorRGB(WHITE));
		s.add_uniform_color("color_scale", get_color_scale());
		if (enable_shadow_maps) {setup_tile_shader_shadow_map(s);}
	}
	if (enable_dlights) {
		setup_dlight_textures(s, 0); // no dlight smap
		s.add_uniform_vector3d("camera_pos", get_camera_pos());
		set_city_lighting_shader_opts(s, lights_bcube, 1, 0); // will reset some values
	}
}


void tile_draw_t::draw_decid_trees(bool reflection_pass, bool shadow_pass) {

	if (to_draw.empty()) return; // nothing to do

	if (to_draw.size() <= 9) { // local shadow computation; almost free, but in most cases this doesn't help
		bool have_trees(0);
		for (auto i = to_draw.begin(); i != to_draw.end() && !have_trees; ++i) {have_trees |= (i->second->num_decid_trees() > 0);}
		if (!have_trees) return;
	}
	bool const enable_billboards(USE_TREE_BILLBOARDS && !shadow_pass);
	bool const enable_shadow_maps(!shadow_pass && shadow_map_enabled()); // && !reflection_pass?
	lod_renderer.resize_zero();
	lod_renderer.set_enabled(enable_billboards); // need full detail rendering in shadow pass, since billboards project poor shadows
	colorRGBA const leaf_color_scale(get_color_scale(0.8, 0.5));

	{ // draw leaves
		bool const leaf_shadow_maps(!(display_mode & 0x0200) && enable_shadow_maps);
		shader_t ls;
		bool const use_fs_smap = 0; // too slow
		if (leaf_shadow_maps) {ls.set_prefix("#define USE_SMAP_SCALE", (use_fs_smap ? 1 : 0));} // VS/FS
		bool const enable_dlights(!reflection_pass && !shadow_pass);
		tree_cont_t::pre_leaf_draw(ls, enable_billboards, shadow_pass, use_fs_smap, leaf_shadow_maps, enable_dlights);
		if (leaf_shadow_maps ) {setup_tile_shader_shadow_map(ls);}
		if (enable_billboards) {lod_renderer.leaf_opacity_loc = ls.get_uniform_loc("opacity");}
		if (!shadow_pass) {set_tree_dither_noise_tex(ls, 1);} // TU=1
		if (!shadow_pass) {ls.add_uniform_color("color_scale", leaf_color_scale);}
		draw_decid_tree_bl(ls, lod_renderer, 0, 1, reflection_pass, shadow_pass, leaf_shadow_maps);
		if (!shadow_pass) {ls.add_uniform_color("color_scale", WHITE);}
		tree_cont_t::post_leaf_draw(ls, shadow_pass);
	}
	{ // draw branches
		shader_t bs;
		bool const enable_dlights(!shadow_pass && !reflection_pass && enable_tree_dlights());
		tree_branch_shader_setup(bs, enable_shadow_maps, 1, shadow_pass, enable_dlights); // enable_opacity=1
		if (!shadow_pass) {set_tree_dither_noise_tex(bs, 1);} // TU=1 (for opacity)
		if (enable_billboards) {lod_renderer.branch_opacity_loc = bs.get_uniform_loc("opacity");}
		draw_decid_tree_bl(bs, lod_renderer, 1, 0, reflection_pass, shadow_pass, enable_shadow_maps);
		bs.add_uniform_vector4d("world_space_offset", vector4d()); // reset
		bs.end_shader();
	}
	lod_renderer.finalize();
	enable_blend(); // for fog transparency

	if (lod_renderer.has_leaves()) { // draw leaf billboards (not in shadow pass)
		shader_t lrs;
		lrs.set_vert_shader("tree_billboard_gs");
		lrs.set_geom_shader("tree_billboard"); // point => 1 quad
		lrs.set_frag_shader("linear_fog.part+ads_lighting.part*+leaf_lighting_comp.part*+noise_dither.part+tree_leaves_billboard");
		billboard_tree_shader_setup(lrs);
		lrs.add_uniform_color("color_scale", leaf_color_scale);
		lrs.set_specular(0.1, 10.0);
		lod_renderer.render_billboards(lrs, 0);
		lrs.clear_specular();
		lrs.end_shader();
	}
	if (lod_renderer.has_branches()) { // draw branch billboards (not in shadow pass)
		shader_t brs;
		brs.set_prefix("#define TREE_BRANCHES", 2); // GS
		brs.set_vert_shader("tree_billboard_gs");
		brs.set_geom_shader("tree_billboard"); // point => 1 quad
		brs.set_frag_shader("linear_fog.part+ads_lighting.part*+noise_dither.part+tree_branches_billboard");
		billboard_tree_shader_setup(brs);
		brs.add_uniform_color("color_scale", get_color_scale());
		brs.add_uniform_vector3d("ref_dir", plus_y);
		lod_renderer.render_billboards(brs, 1);
		brs.end_shader();
	}
	disable_blend();
	if (!reflection_pass) {leaf_color_changed = 0;} // Note: only trees in visible tiles will be updated
}


void tile_draw_t::draw_scenery(bool reflection_pass, bool shadow_pass) {

	bool const enable_shadow_maps(!shadow_pass && shadow_map_enabled());
	shader_t s, vrs;
	vrs.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS - is this needed/useful?

	// draw opaque objects
	// rocks, stumps, and logs are close enough to tree branches that we can reuse the shader
	tree_branch_shader_setup(s, enable_shadow_maps, 0, shadow_pass); // no opacity
	s.begin_shader();
	setup_tt_fog_post(s);
	s.add_uniform_int("tex0", 0);
	for (unsigned i = 0; i < to_draw.size(); ++i) {to_draw[i].second->draw_scenery(s, vrs, 1, 0, reflection_pass, shadow_pass, enable_shadow_maps);}
	tree_scenery_pld.draw_and_clear();
	s.end_shader();

	// draw leaves
	if (enable_shadow_maps) {s.set_prefix("#define USE_SMAP_SCALE", 1);} // FS
	set_leaf_shader(s, 0.9, 0, 0, 1, (shadow_pass ? 0.0 : get_plant_leaf_wind_mag(0)), underwater, enable_shadow_maps, enable_shadow_maps, 1, shadow_pass);
	if (enable_shadow_maps) {setup_tile_shader_shadow_map(s);}
	s.add_uniform_color("color_scale", get_color_scale(1.0, 0.5));
	for (unsigned i = 0; i < to_draw.size(); ++i) {to_draw[i].second->draw_scenery(s, vrs, 0, 1, reflection_pass, shadow_pass, enable_shadow_maps);}
	s.add_uniform_color("color_scale", WHITE);
	s.end_shader();
}

void tile_draw_t::setup_grass_flower_shader(shader_t &s, bool enable_wind, bool use_smap, float dist_const_mult) {

	s.begin_shader();
	if (enable_wind) {setup_wind_for_shader(s, 1);}
	if (use_smap   ) {s.add_uniform_float("smap_atten_cutoff", get_smap_atten_val());}
	s.add_uniform_int("tex0", 0);
	s.add_uniform_int("height_tex", 2);
	s.add_uniform_int("normal_tex", 4);
	s.add_uniform_int("shadow_tex", 6);
	s.add_uniform_float("dist_const", dist_const_mult*get_grass_thresh());
	s.add_uniform_float("dist_slope", 1.0/get_grass_blend_dist());
	setup_tt_fog_post(s);
	setup_cloud_plane_uniforms(s);
	set_tile_xy_vals(s);
	
	if (!grass_exclude1.is_all_zeros()) { // based on primary cube
		s.add_uniform_vector4d("clip_box1", vec4_from_cube_xy(grass_exclude1));
		s.add_uniform_vector4d("clip_box2", vec4_from_cube_xy(grass_exclude2)); // can be zero area
	}
}


// tu's used: 0: grass, 1: wind noise, 2: heightmap, 3: grass weight, 4: normal map, 5: noise, 6: shadow map, 9: cloud noise
void tile_draw_t::draw_grass(bool reflection_pass) {

	if (reflection_pass)    return; // no grass reflections (yet)
	if (player_in_basement) return; // grass can sometimes appear in a building basement, so disable it when the player is in the basement
	//highres_timer_t timer("Draw Grass"); // 0.13/0.38ms
	bool const use_cloud_shadows(GRASS_CLOUD_SHADOWS && cloud_shadows_enabled());
	unsigned num_grass_drawn(0), num_flowers_drawn(0);
	if (use_grass_tess && !check_for_tess_shader()) {use_grass_tess = 0;} // disable tess - not supported

	for (unsigned wpass = 0; wpass < 2; ++wpass) { // wind, no wind
		for (unsigned spass = 0; spass < 2; ++spass) { // shadow maps, no shadow maps
			if (spass == 0 && !shadow_map_enabled()) continue;
			shader_t s;
			bool const enable_wind((display_mode & 0x0100) && wpass == 0);
			bool const enable_tess(use_grass_tess && wpass == 0); // only for nearby grass (use same logic as wind)
			lighting_with_cloud_shadows_setup(s, 0, use_cloud_shadows);
			// Note: when tt_grass_scale_factor is small, we can have a transition from nearby wind to distant within the same tile, so we need to enable height adjust in both cases
			if (wpass == 1 || tt_grass_scale_factor < 1.0) {s.set_prefix("#define DEC_HEIGHT_WHEN_FAR", 0);} // VS
			if (!grass_exclude1.is_all_zeros()) {s.set_prefix("#define ENABLE_VERTEX_CLIP", 0);} // VS - based on primary cube
			//if (!underwater) {s.set_prefix("#define NO_FOG", 1);} // FS - faster, but reduced quality grass/texture blend
			set_smap_enable_for_shader(s, (spass == 0), 0); // VS
			s.set_prefix(make_shader_bool_prefix("enable_grass_wind", enable_wind), 0); // VS
			if (enable_tess) {s.set_prefix("#define NO_FOG_FRAG_COORD", 1);} // FS - needed on some drivers because TC/TE don't have fg_FogFragCoord
			if (!enable_tess) {s.set_prefix("#define NO_GRASS_TESS", 0);} // VS
			s.set_vert_shader("ads_lighting.part*+perlin_clouds.part*+shadow_map.part*+tiled_shadow_map.part*+wind.part*+grass_texture.part+grass_tiled");
			s.set_frag_shader("linear_fog.part+grass_tiled");

			if (enable_tess) {
				s.set_tess_control_shader("grass_tiled");
				s.set_tess_eval_shader("grass_tiled"); // draw calls need to use GL_PATCHES instead of GL_TRIANGLES
				glPatchParameteri(GL_PATCH_VERTICES, 3); // triangles
			}
			setup_grass_flower_shader(s, enable_wind, (spass == 0), 1.0);
			s.add_uniform_int("weight_tex", 3);
			set_noise_tex(s, 5);
			s.add_uniform_float("height", grass_length);
			s.set_specular(0.1, 20.0);
			grass_tile_manager.begin_draw();

			if (enable_tess) {
				s.add_uniform_float("min_tess_level", 1.0);
				s.add_uniform_float("tess_lod_scale", tt_grass_scale_factor);
			}
			int const lt_loc(s.get_attrib_loc("local_translate"));
			enable_instancing_for_shader_loc(lt_loc);

			for (unsigned i = 0; i < to_draw.size(); ++i) {
				tile_t *const tile(to_draw[i].second);
				if (!tile->has_grass()) continue;
				if (tile->using_shadow_maps() != (spass == 0)) continue;
				if ((tile->get_dist_to_camera_in_tiles(0) > 0.5*tt_grass_scale_factor) != (int)wpass) continue; // xyz dist
				num_grass_drawn += tile->draw_grass(s, grass_insts, use_cloud_shadows, enable_tess, lt_loc);
			}
			disable_instancing_for_shader_loc(lt_loc);
			grass_tile_manager.end_draw();
			s.end_shader();
		} // for spass
	} // for wpass

	// draw flowers
	for (unsigned spass = 0; spass < 2; ++spass) { // shadow maps, no shadow maps
		if (flower_density == 0.0) continue; // no flowers
		if (spass == 0 && !shadow_map_enabled()) continue;
		shader_t s;
		lighting_with_cloud_shadows_setup(s, 1, use_cloud_shadows);
		set_smap_enable_for_shader(s, (spass == 0), 1); // FS
		if (!grass_exclude1.is_all_zeros()) {s.set_prefix("#define ENABLE_VERTEX_CLIP", 0);} // VS - based on primary cube
		s.set_vert_shader("texture_gen.part+wind.part*+flowers_tiled");
		//s.set_vert_shader("texture_gen.part+wind.part*+flowers_tiled_gs");
		//s.set_geom_shader("flower_from_pt"); // point => 1 quad
		s.set_frag_shader("linear_fog.part+ads_lighting.part*+perlin_clouds.part*+shadow_map.part*+tiled_shadow_map.part*+flowers_tiled");
		setup_grass_flower_shader(s, 1, (spass == 0), FLOWER_REL_DIST);
		flower_manager_t::setup_flower_shader_post(s);
		enable_blend(); // required for distant flower transition to alpha=0

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			if (to_draw[i].second->using_shadow_maps() != (spass == 0)) continue;
			num_flowers_drawn += to_draw[i].second->draw_flowers(s, use_cloud_shadows);
		}
		disable_blend();
		s.end_shader();
	}
	grass_exclude1.set_to_zeros(); // reset for next frame
	grass_exclude2.set_to_zeros();
	if (DEBUG_TILES) {cout << "grass blades drawn: " << num_grass_drawn << ", flowers drawn: " << num_flowers_drawn << endl;} // up to 2M / 100K
}


void tile_draw_t::draw_animals(bool reflection_pass) {

	shader_t s;
	enable_blend(); // for distance fog

	if (birds_active()) { // day time - applies to birds and butterflies
		// draw birds
		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->draw_birds(s);}
		s.end_shader();

		if (!reflection_pass) { // draw butterflies
			for (unsigned i = 0; i < to_draw.size(); ++i) {to_draw[i].second->draw_bflies(s);}
			s.end_shader();
		}
	}
	if (!reflection_pass) { // draw fish
		for (unsigned i = 0; i < to_draw.size(); ++i) {to_draw[i].second->draw_fish(s);}
		s.end_shader();
	}
	disable_blend(); // for distance fog
}

void tile_draw_t::draw_tile_clouds(bool reflection_pass) { // reflection_pass is unused

	if (!clouds_enabled() || atmosphere < 0.5) return; // only for high atmosphere
	//timer_t timer("Draw Clouds"); // 0.15ms on old computer
	to_draw_clouds.clear();
	unsigned num(0);
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {num += i->second->update_tile_clouds();}
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->get_cloud_draw_list(to_draw_clouds);}
	//cout << "clouds: " << num << endl;
	if (to_draw_clouds.empty()) return;
	sort(to_draw_clouds.begin(), to_draw_clouds.end()); // back-to-front
	vpc_shader_t s; // see draw_scenery()
	//s.set_prefix("#define ENABLE_HEIGHT_ATTEN", 1); // FS Note: not drawn per-tile, and source tile isn't tracked, so not easy to set the correct height texture
	tile_cloud_t::shader_setup(s, 1, 0, -0.12, -0.3, 5, 1, 1); // grayscale, not ridged, with lighting, custom alpha/dist bias, and cloud_mode
	s.add_uniform_float("noise_scale", 0.02);
	s.set_cur_color(WHITE); // unnecessary?
	enable_blend();
	glDepthMask(GL_FALSE); // no depth writing
	set_multisample(0);
	vector3d const xlate(get_tiled_terrain_model_xlate());
	for (auto i = to_draw_clouds.begin(); i != to_draw_clouds.end(); ++i) {i->cloud->draw(s, xlate, i->alpha);}
	set_multisample(1);
	glDepthMask(GL_TRUE);
	disable_blend();
	s.end_shader();
}


void tile_draw_t::update_lightning(bool reflection_pass) {
	if (!reflection_pass) {lightning_strike.update();}
	lightning_strike.draw();
}

void tile_draw_t::clear_vbos_tids() {
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->clear_vbo_tid(nullptr);}
	smap_manager.clear_context();
	free_compute_shader();
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
	return get_tile_from_xy(tile_xy_pair(round_fp(0.5f*(pos.x - (xoff - xoff2)*DX_VAL)/X_SCENE_SIZE), round_fp(0.5f*(pos.y - (yoff - yoff2)*DY_VAL)/Y_SCENE_SIZE)));
}
uint64_t get_tile_id_containing_point(point const &pos) {
	tile_xy_pair const tp(round_fp(0.5f*(pos.x - (xoff - xoff2)*DX_VAL)/X_SCENE_SIZE), round_fp(0.5f*(pos.y - (yoff - yoff2)*DY_VAL)/Y_SCENE_SIZE));
	return (tp.x + (uint64_t(tp.y) << 32));
}
uint64_t get_tile_id_containing_point_no_xyoff(point const &pos) {
	tile_xy_pair const tp(round_fp(0.5*pos.x/X_SCENE_SIZE), round_fp(0.5*pos.y/Y_SCENE_SIZE));
	return (tp.x + (uint64_t(tp.y) << 32));
}

bool tile_draw_t::try_bind_tile_smap_at_point(point const &pos, shader_t &s, bool check_only, unsigned *lod_level) const {
	tile_t const *const tile(get_tile_containing_point(pos));
	return (tile != nullptr && tile->try_bind_shadow_map(s, check_only, lod_level));
}

// Note: region should be less than one tile in size
void tile_draw_t::invalidate_tile_smap_in_region(cube_t region) {
	vector3d const b_ext(get_buildings_max_extent());
	float const road_len_exp(0.5f*get_road_max_len().get_max_val());
	for (unsigned d = 0; d < 2; ++d) {region.expand_in_dim(d, max(road_len_exp, b_ext[d]));} // expand by city tile overlap (should also include trees?)

	for (int y = 0; y < 2; ++y) { // try 4 corners, needed to handle objects that overlap more than one tile
		for (int x = 0; x < 2; ++x) {
			tile_t *const tile(get_tile_containing_point(point(region.d[0][x], region.d[1][y], 0.0)));
			if (tile) {tile->clear_shadow_map(&smap_manager);}
		}
	}
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
bool tile_draw_t::check_cube_int_trees(cube_t const &c) const {
	tile_t const *const tile(get_tile_containing_point(c.get_cube_center())); // assumes cube is contained in one tile
	return (tile ? tile->check_cube_int_trees(c) : 0);
}

int tile_draw_t::get_tid_under_point(point const &pos) const {
	tile_t const *const tile(get_tile_containing_point(pos));
	return (tile ? tile->get_tid_under_point(pos) : -1);
}


bool tile_draw_t::line_intersect_mesh(point const &v1, point const &v2, float &t, tile_t *&intersected_tile, int &xpos, int &ypos, float inc_trees) const {

	t = 2.0; // > 1.0
	intersected_tile = nullptr;

	// Note: could start with get_tile_containing_point(v1) and walk the line to v2 using DDA, etc.
	for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
		float tn(1.0);
		int new_xpos(0), new_ypos(0);

		// Note: could make this faster by passing tmin, tmax into line_intersect_mesh() here and using for early termination,
		// but this code is plenty fast enough to do a single query each frame as it is
		if (i->second->line_intersect_mesh(v1, v2, tn, new_xpos, new_ypos, inc_trees) && tn < t) {
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
	// Note: it would be better if we could just translate the shadow map when the scene shifts, but this seems fairly complex to track and get right
	int const new_dxoff(xoff - xoff2), new_dyoff(yoff - yoff2);
	bool const new_off(new_dxoff != dxoff || new_dyoff != dyoff);
	dxoff = new_dxoff; dyoff = new_dyoff;
	return (smap_data_t::needs_update(lpos) || new_off);
}


tile_t *get_tile_from_xy  (tile_xy_pair const &tp) {return terrain_tile_draw.get_tile_from_xy(tp);}
float update_tiled_terrain(float &min_camera_dist) {return terrain_tile_draw.update(min_camera_dist);}
void pre_draw_tiled_terrain() {terrain_tile_draw.pre_draw();}
void show_tiled_terrain_debug_stats() {terrain_tile_draw.show_debug_stats(0);} // calc_mem_only=0
uint64_t get_tiled_terrain_gpu_mem() {return terrain_tile_draw.show_debug_stats(1);} // calc_mem_only=1


colorRGBA get_inf_terrain_mod_color() {

	colorRGBA const colors[NUM_FIRE_MODES] = {WHITE, GREEN, RED, BLUE, YELLOW, CYAN, BROWN, WHITE};
	assert(inf_terrain_fire_mode < NUM_FIRE_MODES);
	return colors[inf_terrain_fire_mode];
}

bool line_intersect_tiled_mesh_get_tile(point const &v1, point const &v2, point &p_int, tile_t *&tile, bool inc_trees) {

	float t(0.0);
	int xpos(0), ypos(0); // unused
	if (!terrain_tile_draw.line_intersect_mesh(v1, v2, t, tile, xpos, ypos, inc_trees)) return 0;
	p_int = v1 + t*(v2 - v1);
	return 1;
}

int get_tiled_terrain_tid_under_point(point const &pos) {
	return terrain_tile_draw.get_tid_under_point(pos);
}

void draw_brush_shape(float xval, float yval, float radius, float z1, float z2, bool is_square) {

	if (is_square) { // square, projected to cube
		draw_cube(point(xval, yval, 0.5f*(z1 + z2)), 2.0f*radius, 2.0f*radius, z2-z1, 0);
	}
	else { // circle, projected to cylinder
		draw_cylinder_at(point(xval, yval, z1), z2-z1, radius, radius, N_CYL_SIDES, 1); // with ends
	}
}

void render_tt_models(int reflection_pass, bool transparent_pass) {

	if (reflection_pass && !enable_tt_model_reflect) return;
	vector3d const xlate(get_tiled_terrain_model_xlate());
	render_models(0, reflection_pass, (transparent_pass ? 2 : 1), xlate);
}

void draw_tiled_terrain(int reflection_pass) {

	//RESET_TIME;
	// skip terrain draw if the player is in a closed elevator so that we don't see the terrain when crossing the ground;
	// but the terrain may be visible through tw window through the crack between the doors, though this isn't too much of an issue
	if (player_in_elevator >= 2) return;
	bool const disable_dclamp(enable_depth_clamp &&  camera_in_building && !reflection_pass   ); // helps with terrain covering basement stairs entrance
	bool const enable_dclamp(!enable_depth_clamp && !camera_in_building && camera_surf_collide); // helps prevent camera clipping through the terrain
	if (disable_dclamp) {glDisable(GL_DEPTH_CLAMP);}
	if (enable_dclamp ) {glEnable (GL_DEPTH_CLAMP);}
	// don't need to draw the bottom of the terrain when in the basement; for some reason the faces are backwards
	if (player_in_basement)  {glEnable (GL_CULL_FACE); glCullFace(GL_FRONT);}
	terrain_tile_draw.draw(reflection_pass);
	if (player_in_basement)  {glDisable(GL_CULL_FACE); glCullFace(GL_BACK);}
	if (disable_dclamp) {glEnable (GL_DEPTH_CLAMP);} // restore
	if (enable_dclamp)  {glDisable(GL_DEPTH_CLAMP);} // restore
	//glFinish(); PRINT_TIME("Tiled Terrain Draw"); //exit(0);
	if (reflection_pass) return; // nothing else to do

	if (inf_terrain_fire_mode != FM_NONE) { // use a bool instead?
		point const v1(get_camera_pos()), v2(v1 + cview_dir*FAR_CLIP);
		point hit_pos;
		tile_t *tile(nullptr);

		if (line_intersect_tiled_mesh_get_tile(v1, v2, hit_pos, tile, 0)) { // inc_trees=0
			draw_single_colored_sphere(hit_pos, 0.1, N_SPHERE_DIV, RED);

			// modification marker area rendered with a cylinder + stencil test to mask to mesh surface
			assert(tile != nullptr);
			float const tzmin(min(hit_pos.z, tile->get_zmin())), tzmax(max(hit_pos.z, tile->get_tile_zmax())), dz(tzmax - tzmin);
			float const z1(tzmin - dz - SMALL_NUMBER), z2(tzmax + dz + SMALL_NUMBER);
			float const radius((cur_brush_param.get_radius() + 0.5)*HALF_DXY); // do we need a nonuniform x/y scale for non-square meshes?
			bool const is_square(cur_brush_param.shape == BSHAPE_CONST_SQ);
			enable_blend();
			shader_t s;
			s.begin_color_only_shader(colorRGBA(get_inf_terrain_mod_color(), 0.2)); // make mostly transparent (wireframe?)
			setup_stencil_buffer_write();
			glStencilOpSeparate(GL_FRONT, GL_INCR_WRAP, GL_INCR_WRAP, GL_KEEP);
			glStencilOpSeparate(GL_BACK,  GL_DECR_WRAP, GL_DECR_WRAP, GL_KEEP);
			draw_brush_shape(hit_pos.x, hit_pos.y, radius, z1, z2, is_square);
			glStencilFunc(GL_NOTEQUAL, -1, ~0U); // one front face only
			glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			draw_brush_shape(hit_pos.x, hit_pos.y, radius, z1, z2, is_square);
			end_stencil_write();
			disable_blend();
			s.end_shader();
		}
	}
	if (camera_surf_collide) {
		//int const tid(get_tiled_terrain_tid_under_point(get_camera_pos())); // TESTING
	}
}

void draw_tiled_terrain_lightning(bool reflection_pass) {terrain_tile_draw.update_lightning(reflection_pass);}
void end_tiled_terrain_lightning() {terrain_tile_draw.end_lightning();}
void clear_tiled_terrain(bool no_regen_buildings) {terrain_tile_draw.clear(no_regen_buildings);}
void draw_tiled_terrain_clouds(bool reflection_pass) {terrain_tile_draw.draw_tile_clouds(reflection_pass);}
void draw_tiled_terrain_decid_tree_shadows() {terrain_tile_draw.draw_decid_tree_shadows();}
void reset_tiled_terrain_state() {terrain_tile_draw.clear_vbos_tids();}
void clear_tiled_terrain_shaders() {terrain_tile_draw.free_compute_shader();}
void draw_tiled_terrain_water(shader_t &s, float zval) {terrain_tile_draw.draw_water(s, zval);}
bool check_player_tiled_terrain_collision() {return terrain_tile_draw.check_player_collision();}
bool sphere_int_tiled_terrain(point &pos, float radius) {return terrain_tile_draw.check_sphere_collision(pos, radius);}
bool cube_int_tiled_terrain_trees(cube_t const &c) {return terrain_tile_draw.check_cube_int_trees(c);}
float get_tiled_terrain_water_level() {return (is_water_enabled() ? water_plane_z : terrain_tile_draw.get_actual_zmin());}

bool try_bind_tile_smap_at_point(point const &pos, shader_t &s, bool check_only, unsigned *lod_level) {
	return terrain_tile_draw.try_bind_tile_smap_at_point(pos, s, check_only, lod_level);
}
// defer update until tile draw (if called from non-drawing thread); region and pos are in camera space
void invalidate_tile_smap_in_region(cube_t const &region, bool repeat_next_frame) {tile_smaps_to_clear.emplace_back(region, repeat_next_frame);}

void invalidate_tile_smap_at_pt(point const &pos, float radius, bool repeat_next_frame) {
	cube_t region;
	region.set_from_sphere(pos, radius);
	invalidate_tile_smap_in_region(region, repeat_next_frame);
}

// *** tree/grass addition/removal ***


void tile_draw_t::add_or_remove_trees_at(point const &pos, float radius, bool add_trees, int brush_shape) {

	static point last_pos(all_zeros);
	static float last_radius(0.0);

	if (add_trees) {
		if (pos == last_pos && radius == last_radius) return; // brush applied at same point (mouse not moving/zero delay), ignore it (weight is unused)
		last_pos = pos; last_radius = radius;
	}
	//timer_t timer("Add/Remove Trees"); // add pine=4.4 add decid=3.5
	cube_t update_bcube(all_zeros);
	vector<tile_t *> near_tiles;

	// Note: not suitable for openmp because it modifies shared state (smap_manager, near_tiles, shared_tree_data, VBOs, xoff2, yoff2)
	// also, many edit operations will affect a single tile anyway, which won't distribute well; and this step is only part of the CPU time
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
		if (i->second->add_or_remove_trees_at(pos, radius, add_trees, brush_shape, smap_manager, update_bcube)) {near_tiles.push_back(i->second.get());}
	}
	if (update_bcube.is_all_zeros()) return; // no trees updated

	for (auto i = near_tiles.begin(); i != near_tiles.end(); ++i) { // update shadows for nearby tiles
		if ((*i)->get_mesh_bcube().intersects(update_bcube)) {(*i)->register_tree_change(smap_manager);}
	}
}

void tile_draw_t::add_or_remove_grass_at(point const &pos, float radius, bool add_grass, int brush_shape, float brush_weight) {
	//timer_t timer("Add/Remove Grass"); // rem=1.4 add=1.4
	for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {i->second->add_or_remove_grass_at(pos, radius, add_grass, brush_shape, brush_weight);}
}

void update_trees_bcube(point const &tpos, float tradius, cube_t &bcube) {
	if (bcube.is_all_zeros()) {bcube.set_from_sphere(tpos, tradius);} else {bcube.union_with_sphere(tpos, tradius);}
}
template <typename T> bool remove_tree(vector<T> &v, unsigned &i, point const &pos, float rradius, bool is_square, point const &xlate, cube_t &rem_bcube) {

	point const &tpos(v[i].get_center());
	if (fabs(tpos.x - pos.x) > rradius || fabs(tpos.y - pos.y) > rradius) return 0;
	if (!is_square && !dist_xy_less_than(tpos, pos, rradius)) return 0;
	update_trees_bcube(tpos+xlate, 2.0*v[i].get_radius(), rem_bcube);
	remove_element(v, i);
	return 1;
}

void tile_t::register_tree_change(tile_shadow_map_manager &smap_manager) {
	clear_shadow_map(&smap_manager);
	tree_map.clear(); // regenerate tree AO shadows
	clear_shadows();
	recalc_tree_grass_weights = 1; // force recreation of tile texture, even if there's no grass
}

bool tile_t::mesh_sphere_intersect(point const &pos, float rradius) const {
	if (!dist_less_than(pos, get_center(), (radius + rradius))) return 0; // tile not within sphere's radius
	return sphere_cube_intersect(pos, rradius, get_mesh_bcube()); // tighter intersection test
}

bool tile_t::use_as_occluder() const {
	return (!has_tunnel && rel_dist_to_camera_xy_lt(OCCLUDER_DIST) && is_visible());
}

template <typename T> bool tile_t::add_new_trees(T &trees, tile_offset_t const &toff, cube_t &update_bcube, float &tzmax, point const &tpos, float rradius, bool is_square) {

	unsigned const start_sz(trees.size());
	int const orig_xoff2(xoff2), orig_yoff2(yoff2);
	xoff2 = -toff.dxoff; yoff2 = -toff.dyoff; // translate so that trees are generated relative to toff, rather than current offset
	trees.gen_trees_tt_within_radius(x1+toff.dxoff, y1+toff.dyoff, x2+toff.dxoff, y2+toff.dyoff, tpos, rradius, is_square, mesh_dz, this);
	xoff2 = orig_xoff2; yoff2 = orig_yoff2; // translate back
	postproc_trees(trees, tzmax);
	vector3d const xlate(toff.get_xlate());

	for (unsigned i = start_sz; i < trees.size(); ++i) {
		update_trees_bcube(trees[i].get_center()+xlate, 2.0*trees[i].get_radius(), update_bcube);
	}
	return (trees.size() > start_sz);
}

// return value: 0 = not within sphere, 1 = no trees updated, 2 = trees updated
int tile_t::add_or_remove_trees_at(point const &pos, float rradius, bool add_trees, int brush_shape, tile_shadow_map_manager &smap_manager, cube_t &update_bcube) {
	
	if (!mesh_sphere_intersect(pos, (1.1*rradius + 2.0*trmax))) return 0; // sphere not near this tile, and not overlapping any tile trees
	if (!mesh_sphere_intersect(pos, rradius)) return 1; // not overlapping this tile, but close
	bool pine_changed(0), decid_changed(0);
	vector3d const pine_xlate(ptree_off.get_xlate()), decid_xlate(dtree_off.get_xlate());
	point const pt_pos(pos - pine_xlate), dt_pos(pos - decid_xlate);
	bool const is_square(brush_shape == BSHAPE_CONST_SQ);

	// remove trees first, before adding, so that there are no overlaps of new trees and old trees
	for (unsigned i = 0; i < pine_trees.size();  ++i) {pine_changed  |= remove_tree(pine_trees,  i, pt_pos, rradius, is_square, pine_xlate,  update_bcube);}
	for (unsigned i = 0; i < decid_trees.size(); ++i) {decid_changed |= remove_tree(decid_trees, i, dt_pos, rradius, is_square, decid_xlate, update_bcube);}

	if (add_trees) { // what if trees are not generated? can we get here in that case?
		if (can_have_pine_palm_trees() && pine_trees_generated ()) {pine_changed  |= add_new_trees(pine_trees,  ptree_off, update_bcube, ptzmax, pt_pos, rradius, is_square);}
		if (can_have_decid_trees() && decid_trees.was_generated()) {decid_changed |= add_new_trees(decid_trees, dtree_off, update_bcube, dtzmax, dt_pos, rradius, is_square);}
	}
	if (!pine_changed && !decid_changed) return 1; // no trees updated
	if (pine_changed) {clear_pine_tree_vbos();}
	register_tree_change(smap_manager);
	return 2; // trees updated
}

bool tile_t::add_or_remove_grass_at(point const &pos, float rradius, bool add_grass, int brush_shape, float brush_weight) {

	if (rradius == 0.0) return 0;
	if (weight_data.empty()) return 0; // texture weights not yet generated, can this happen?
	if (!mesh_sphere_intersect(pos, rradius)) return 0; // not overlapping this tile
	//if (!add_grass && gen_grass_map() && grass_blocks.empty()) return 0; // known to have no grass (is this safe? does it help?)
	bool updated(0);
	unsigned const tsize(stride), num_texels(tsize*tsize);
	assert(weight_data.size() == 4*num_texels);
	int sand_tex_ix(-1), dirt_tex_ix(-1), grass_tex_ix(-1), rock_tex_ix(-1), snow_tex_ix(-1);
	get_texture_ixs(sand_tex_ix, dirt_tex_ix, grass_tex_ix, rock_tex_ix, snow_tex_ix);
	bool const is_square(brush_shape == BSHAPE_CONST_SQ);
	unsigned const grass_block_dim(get_grass_block_dim());
	float const r_inv(1.0/rradius), dz_inv(1.0f/(zmax - zmin)), xy_mult(1.0/float(size));
	float const bweight(10.0*brush_weight); // normal brush weight is 0.001 to 0.512
	float const llc_x(get_xval(x1 + xoff - xoff2)), llcy(get_yval(y1 + yoff - yoff2));
	point pt(llc_x, llcy, 0.0); // z is unused
	unsigned xl(size), yl(size), xh(0), yh(0); // for add_grass only; start as denormalized range

	for (unsigned y = 0; y < tsize; ++y, pt.y += DY_VAL) {
		if (fabs(pt.y - pos.y) > rradius) continue;
		pt.x = llc_x;

		for (unsigned x = 0; x < tsize; ++x, pt.x += DX_VAL) {
			if (fabs(pt.x - pos.x) > rradius) continue;
			if (!is_square && !dist_xy_less_than(pt, pos, rradius)) continue; // check for round shapes
			unsigned const off(4*(y*tsize + x));
			unsigned char &grass_weight(weight_data[off + grass_tex_ix]);
			if (grass_weight == (add_grass ? 255 : 0)) continue; // full/no grass here
			float delta(bweight);
			adjust_brush_weight(delta, p2p_dist_xy(pt, pos)*r_inv, brush_shape);
			unsigned const ix(y*zvsize + x);
			float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
			float const mhmin(min(min(mh00, mh01), min(mh10, mh11))), mhmax(max(max(mh00, mh01), max(mh10, mh11)));

			if (add_grass) {
				if (delta >= 0.99) { // full addition
					grass_weight = 255; // max grass
					UNROLL_4X(if (i_ != grass_tex_ix) {weight_data[off+i_] = 0;}); // set all other weights to 0
				}
				else if (delta > 0.01) { // partial addition of grass
					unsigned char const prev_gw(grass_weight);
					grass_weight = (unsigned char)min(255.0f, (grass_weight + 255.0f*delta));
					unsigned char grass_added(grass_weight - prev_gw);

					for (unsigned i = 0; i < 4 && grass_added > 0; ++i) { // remove other weights until weights are balanced
						if ((int)i == grass_tex_ix) continue; // skip grass
						if (weight_data[off+i] == 0) continue; // none of this layer
						unsigned char const num_rem(min(weight_data[off+i], grass_added)); // should depend on ratio of other weights (but rounding is difficult)
						weight_data[off+i] -= num_rem;
						grass_added        -= num_rem;
					}
				}
				add_grass_block_at(x, y, mhmin, mhmax, grass_block_dim);
				xl = min(xl, x); xh = max(xh, min(x+1, size));
				yl = min(yl, y); yh = max(yh, min(y+1, size));
				updated = 1;
			}
			else { // remove grass
				unsigned char const prev_gw(grass_weight);
				grass_weight = (unsigned char)max(0.0f, (grass_weight - 255.0f*delta));
				unsigned char grass_rem(prev_gw - grass_weight);
				if (grass_rem == 0) continue; // no change
				int k1(0), k2(0);
				float t(0.0);
				get_tids((relh_adj_tex + (mhmax - zmin)*dz_inv), k1, k2, &t);
				unsigned char grass_rem2(t*grass_rem), grass_rem1(grass_rem - grass_rem2);
				if (k1 == grass_tex_ix) {k1 = dirt_tex_ix;} // replace grass with dirt
				if (k2 == grass_tex_ix) {k2 = dirt_tex_ix;} // replace grass with dirt
				if (k2 < 4) {weight_data[off + k2] += grass_rem2;} // not snow
				if (k1 < 4) {weight_data[off + k1] += grass_rem1;} // not snow
				float const dirt_scale(BILINEAR_INTERP(params, dirt, float(x)*xy_mult, float(y)*xy_mult));

				if (dirt_scale < 1.0) { // apply dirt scale: convert dirt to sand
					unsigned char &dirt_w(weight_data[off + dirt_tex_ix]);
					weight_data[off + sand_tex_ix] += (unsigned char)((1.0 - dirt_scale)*dirt_w);
					dirt_w = (unsigned char)(dirt_scale*dirt_w);
				}
				updated = 1;
			}
		} // for x
	} // for y
	if (!updated) return 0;

	if (add_grass) {
		//flowers.clear(); // clear and regenerate flowers (simple but slow)
		flowers.update_subrange(weight_data, stride, x1-xoff2, y1-yoff2, xl, yl, xh, yh); // clear and regenerate range (faster)
	}
	else if (!flowers.empty()) { // remove flowers under grass if grass has been removed
		point const flower_xlate(get_xval(x1 + xoff - xoff2), get_yval(y1 + yoff - yoff2), 0.0);
		flowers.clear_within((pos - flower_xlate), rradius, is_square);
	}
	if (!add_grass && has_grass()) { // increases edit time slightly but decreased draw time slightly
		bool has_grass(0);

		for (unsigned i = 0; i < num_texels; ++i) {
			if (weight_data[4*i + grass_tex_ix] > 0) {has_grass = 1; break;}
		}
		if (!has_grass) {grass_blocks.clear();} // clear all grass blocks when there is no more grass
	}
	create_or_update_weight_tex();
	calc_avg_mesh_color(); // optional
	return 1;
}


// *** heightmap modification and queries ***


bool line_intersect_tiled_mesh(point const &v1, point const &v2, point &p_int, bool inc_trees) {
	tile_t *tile(nullptr); // unused
	return line_intersect_tiled_mesh_get_tile(v1, v2, p_int, tile, inc_trees);
}

void change_inf_terrain_fire_mode(int val, bool mouse_wheel) {

	if (have_buildings()) { // no terrain editing when there are buildings as this doesn't work properly
		building_gameplay_action_key(((val > 0) ? 1 : 0), mouse_wheel); // use this for gameplay instead
		return;
	}
	inf_terrain_fire_mode = (inf_terrain_fire_mode + NUM_FIRE_MODES + val) % NUM_FIRE_MODES;
	
	if (!using_tiled_terrain_hmap_tex() && inf_terrain_fire_mode >= FM_INC_MESH && inf_terrain_fire_mode <= FM_FLATTEN) {
		inf_terrain_fire_mode = ((val > 0) ? (unsigned)FM_REM_TREES : (unsigned)FM_NONE); // skip over mesh heightmap edit modes, which won't work in this case
	}
	string const modes[NUM_FIRE_MODES] = {"Look Only", "Increase Mesh Height", "Decrease Mesh Height", "Flatten Mesh", "Remove Trees", "Add Trees", "Remove Grass", "Add Grass"};
	print_text_onscreen(modes[inf_terrain_fire_mode], WHITE, 1.0, TICKS_PER_SECOND, 1); // 1 second
	play_switch_weapon_sound();
}

tile_t *get_tile_for_xy(int x, int y) {

	int const tsz(get_tile_size());
	if (x < 0) {x -= tsz-1;} // handle truncation toward lower integer
	if (y < 0) {y -= tsz-1;}
	return get_tile_from_xy(tile_xy_pair(x/tsz, y/tsz));
}

void inf_terrain_fire_weapon() {

	if (inf_terrain_fire_mode == FM_NONE) {
		tt_fire_button_down = 1;
		return; // ignore
	}
	//RESET_TIME;
	static double last_tfticks(0.0);
	if ((tfticks - last_tfticks) <= cur_brush_param.delay) return; // limit firing rate
	last_tfticks = tfticks;
	point const v1(get_camera_pos()), v2(v1 + cview_dir*FAR_CLIP);
	unsigned const bradius(cur_brush_param.get_radius());
	point p_int;
	
	// FM_NONE, FM_INC_MESH, FM_DEC_MESH, FM_FLATTEN, FM_REM_TREES, FM_ADD_TREES, FM_REM_GRASS, FM_ADD_GRASS
	if (inf_terrain_fire_mode == FM_REM_TREES || inf_terrain_fire_mode == FM_ADD_TREES) { // tree addition/removal
		if (line_intersect_tiled_mesh(v1, v2, p_int)) {
			terrain_tile_draw.add_or_remove_trees_at(p_int, (bradius + 0.5)*HALF_DXY, (inf_terrain_fire_mode == FM_ADD_TREES), cur_brush_param.shape);
		}
		return;
	}
	if (inf_terrain_fire_mode == FM_REM_GRASS || inf_terrain_fire_mode == FM_ADD_GRASS) { // grass addition/removal
		if (grass_density == 0) return; // disabled

		if (line_intersect_tiled_mesh(v1, v2, p_int)) {
			terrain_tile_draw.add_or_remove_grass_at(p_int, (bradius + 0.5)*HALF_DXY, (inf_terrain_fire_mode == FM_ADD_GRASS), cur_brush_param.shape, cur_brush_param.get_delta_mag());
		}
		return;
	}
	float t(0.0); // unused
	tile_t *tile(NULL);
	int xpos(0), ypos(0);
	if (!terrain_tile_draw.line_intersect_mesh(v1, v2, t, tile, xpos, ypos, 0)) return; // inc_trees=0
	// Note: update is slow when trees are enabled
	unsigned shape(cur_brush_param.shape);
	if (inf_terrain_fire_mode == FM_FLATTEN) {shape = ((shape == BSHAPE_CONST_SQ) ? (unsigned)BSHAPE_FLAT_SQ : (unsigned)BSHAPE_FLAT_CIR);} // enable a flattening shape
	float const delta_mag(cur_brush_param.get_delta_mag()*((inf_terrain_fire_mode == FM_INC_MESH) ? 1.0 : -1.0));
	tex_mod_map_manager_t::hmap_val_t const base_delta(terrain_hmap_manager.scale_delta(delta_mag));
	tex_mod_map_manager_t::hmap_brush_t const brush(xpos, ypos, base_delta, bradius, shape);
	terrain_hmap_manager.apply_brush(brush, tile, 1); // cache
	//PRINT_TIME("Hmap Brush");
}

void inf_terrain_undo_hmap_mod() {

	if (!inf_terrain_fire_mode) return; // ignore
	// Note: doesn't handle tree and grass edits, only terrain height edits
	if (!using_tiled_terrain_hmap_tex()) return; // height editing only supported for hmap terrain
	tex_mod_map_manager_t::hmap_brush_t brush;
	if (!terrain_hmap_manager.pop_last_brush(brush)) return;
	if (brush.is_flatten_brush()) return; // can't undo this brush since it's lossy
	// Note: won't work if clamping to min/max height occurred when applying the brush the first time
	brush.delta = -brush.delta; // invert
	terrain_hmap_manager.apply_brush(brush, get_tile_for_xy(brush.x, brush.y), 0); // don't cache
}

void flatten_hmap_region(cube_t const &cube) {
	if (using_tiled_terrain_hmap_tex()) {terrain_hmap_manager.flatten_region(cube);}
}

void write_heightmap_png(string const &fn) {terrain_hmap_manager.write_png(fn);}


