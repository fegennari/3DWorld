// 3D World - Main function, glut callbacks, mouse and keyboard processing, window management, config file variable handling, etc.
// by Frank Gennari
// 3/10/02

#include "mesh.h"
#include "main.h"
#include "sinf.h"
#include "ship.h"
#include "ship_util.h"
#include "player_state.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "gl_ext_arb.h"
#include "model3d.h"
#include "openal_wrap.h"
#include "file_utils.h"
#include "draw_utils.h"
#include "tree_leaf.h"
#include <set>
#include <thread> // for std::thread::hardware_concurrency()

#ifdef _WIN32 // wglew.h seems to be Windows only
#include <GL/wglew.h> // for wglSwapIntervalEXT
#else // linux
#include <GL/glxew.h>
#endif

#ifdef _WIN32
extern "C" { // force use of dedicated GPUs rather than Intel integrated graphics
	__declspec(dllexport) DWORD NvOptimusEnablement = 1;
	__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

using namespace std;
typedef set<unsigned char>::iterator keyset_it;


int const INIT_DMODE       = 0x010F;
bool const MOUSE_LOOK_DEF  = 0;
int const STARTING_INIT_X  = 0; // setting to 1 seems safer but less efficient
float const DEF_CRADIUS    = 10.0;
float const DEF_CTHETA     = -1.0;
float const DEF_CPHI       = 1.5;
float const DEF_UPTHETA    = 0.0;
float const DEF_CAMY       = 1.0;
float const MAP_SHIFT      = 20.0;
float const MAP_ZOOM       = 1.5;
float const D_TIMESTEP     = 1.5;
float const MOUSE_R_ADJ    = 0.01;   // mouse zoom radius adjustment
float const MOUSE_TRAN_ADJ = 0.01;   // mouse translate adjustment
float const MOUSE_ANG_ADJ  = PI/350; // mouse camera angle increment
float const MA_TOLERANCE   = 0.0001; // tolerance adjustment
float const CAMERA_AIR_CONT= 0.5;
float const DEF_OCEAN_WAVE_HEIGHT = 0.01;
float const DEF_CAMERA_RADIUS     = 0.06;
int const START_MODE       = WMODE_GROUND; // 0 = standard mesh, 1 = planet/universe, 2 = infinite terrain


char const *const defaults_file  = "defaults.txt";
char const *const dstate_file    = "state.txt";
char const *const dmesh_file     = "mesh.txt";
char const *const dcoll_obj_file = "coll_objs/coll_objs.txt";
char const *const dship_def_file = "universe/ship_defs.txt";
char *state_file(nullptr), *mesh_file(nullptr), *coll_obj_file(nullptr);
char *mh_filename(nullptr), *mh_filename_tt(nullptr), *mesh_diffuse_tex_fn(nullptr), *ship_def_file(nullptr), *snow_file(nullptr);
char *lighting_file[NUM_LIGHTING_TYPES] = {0};


// Global Variables - most of these are set by the config file reader and used in other files.
// I decided to use global variables here rather than a global config class to avoid frequent recompile of all code
// every time a config option is added/changed, because almost every file would need to include the class definition/header.
// Note that these are all the default values when no config variable is specified.
bool combined_gu(0), underwater(0), kbd_text_mode(0), univ_stencil_shadows(1), use_waypoint_app_spots(0), enable_tiled_mesh_ao(0), tiled_terrain_only(0);
bool show_lightning(0), disable_shader_effects(0), use_waypoints(0), group_back_face_cull(0), start_maximized(0), claim_planet(0), skip_light_vis_test(0);
bool no_smoke_over_mesh(0), enable_model3d_tex_comp(0), global_lighting_update(0), lighting_update_offline(0), mesh_difuse_tex_comp(1), smoke_dlights(0), keep_keycards_on_death(0);
bool texture_alpha_in_red_comp(0), use_model3d_tex_mipmaps(1), mt_cobj_tree_build(0), two_sided_lighting(0), inf_terrain_scenery(1), invert_model_nmap_bscale(0);
bool gen_tree_roots(1), fast_water_reflect(0), vsync_enabled(0), use_voxel_cobjs(0), disable_sound(0), enable_depth_clamp(0), volume_lighting(0), no_subdiv_model(0);
bool detail_normal_map(0), init_core_context(0), use_core_context(0), enable_multisample(1), dynamic_smap_bias(0), model3d_wn_normal(0), snow_shadows(0), user_action_key(0);
bool enable_dlight_shadows(1), tree_indir_lighting(0), ctrl_key_pressed(0), only_pine_palm_trees(0), enable_gamma_correct(0), use_z_prepass(0), reflect_dodgeballs(0);
bool store_cobj_accum_lighting_as_blocked(0), all_model3d_ref_update(0), begin_motion(0), enable_mouse_look(MOUSE_LOOK_DEF), enable_init_shields(1), tt_triplanar_tex(0);
bool enable_model3d_bump_maps(1), use_obj_file_bump_grayscale(1), invert_bump_maps(0), use_interior_cube_map_refl(0), enable_cube_map_bump_maps(1), no_store_model_textures_in_memory(0);
bool enable_model3d_custom_mipmaps(1), flatten_tt_mesh_under_models(0), smileys_chase_player(0), disable_fire_delay(0), disable_recoil(0), mesh_size_locked(0);
bool enable_dpart_shadows(0), enable_tt_model_reflect(1), enable_tt_model_indir(0), auto_calc_tt_model_zvals(0), use_model_lod_blocks(0), enable_translocator(0), enable_grass_fire(0);
bool disable_model_textures(0), start_in_inf_terrain(0), allow_shader_invariants(1), config_unlimited_weapons(0), disable_tt_water_reflect(0), allow_model3d_quads(1);
bool enable_timing_profiler(0), fast_transparent_spheres(0), force_ref_cmap_update(0), use_instanced_pine_trees(0), enable_postproc_recolor(0), draw_building_interiors(0);
bool toggle_room_light(0), teleport_to_screenshot(0), merge_model_objects(0), reverse_3ds_vert_winding_order(1), disable_dlights(0), voxel_add_remove(0), enable_ground_csm(0);
bool enable_hcopter_shadows(0), pre_load_full_tiled_terrain(0), disable_blood(0), enable_model_animations(1), rotate_trees(0), invert_model3d_faces(0), play_gameplay_alert(1);
bool player_custom_start_pos(0);
int xoff(0), yoff(0), xoff2(0), yoff2(0), rand_gen_index(0), mesh_rgen_index(0), camera_change(1), camera_in_air(0), auto_time_adv(0);
int animate(1), animate2(1), draw_model(0), init_x(STARTING_INIT_X), fire_key(0), do_run(0), init_num_balls(-1), change_wmode_frame(0);
int game_mode(0), map_mode(0), load_hmv(0), load_coll_objs(1), read_landscape(0), screen_reset(0), mesh_seed(0), rgen_seed(1);
int display_framerate(1), init_resize(1), temp_change(0), is_cloudy(0), recreated(1), cloud_model(0), force_tree_class(-1);
int invert_mh_image(0), voxel_editing(0), min_time(0), show_framerate(0), preproc_cube_cobjs(0), use_voxel_rocks(2);
int camera_view(0), camera_reset(1), camera_mode(0), camera_surf_collide(1), camera_coll_smooth(0), use_smoke_for_fog(0);
int window_width(512), window_height(512), map_color(1); // window dimensions, etc.
int border_height(20), border_width(4), world_mode(START_MODE), display_mode(INIT_DMODE), do_read_mesh(0);
int last_mouse_x(0), last_mouse_y(0), m_button(0), mouse_state(1), maximized(0), fullscreen(0), verbose_mode(0), leaf_color_changed(0);
int do_zoom(0), disable_universe(0), disable_inf_terrain(0), precip_mode(0), building_action_key(0);
int num_trees(0), num_smileys(1), srand_param(3), left_handed(0), mesh_scale_change(0);
int pause_frame(0), show_fog(0), spectate(0), b2down(0), free_for_all(0), teams(2), show_scores(0), universe_only(0);
int reset_timing(0), read_heightmap(0), default_ground_tex(-1), num_dodgeballs(1), INIT_DISABLE_WATER, ground_effects_level(2);
int enable_fsource(0), run_forward(0), advanced(0), dynamic_mesh_scroll(0), default_anim_id(-1), stats_display_mode(0);
int read_snow_file(0), write_snow_file(0), mesh_detail_tex(NOISE_TEX);
int read_light_files[NUM_LIGHTING_TYPES] = {0}, write_light_files[NUM_LIGHTING_TYPES] = {0};
unsigned num_snowflakes(0), create_voxel_landscape(0), hmap_filter_width(0), num_dynam_parts(100), snow_coverage_resolution(2), show_map_view_fractal(0);
unsigned num_birds_per_tile(2), num_fish_per_tile(15), num_bflies_per_tile(4);
unsigned erosion_iters(0), erosion_iters_tt(0), skybox_tid(0), tiled_terrain_gen_heightmap_sz(0), game_mode_disable_mask(0), num_frame_draw_calls(0);
float NEAR_CLIP(DEF_NEAR_CLIP), FAR_CLIP(DEF_FAR_CLIP), system_max_orbit(1.0), sky_occlude_scale(0.0), tree_slope_thresh(5.0), mouse_sensitivity(1.0), tt_grass_scale_factor(1.0);
float water_plane_z(0.0), base_gravity(1.0), crater_depth(1.0), crater_radius(1.0), disabled_mesh_z(FAR_CLIP), vegetation(1.0), atmosphere(1.0), biome_x_offset(0.0);
float mesh_file_scale(1.0), mesh_file_tz(0.0), speed_mult(1.0), mesh_z_cutoff(-FAR_CLIP), relh_adj_tex(0.0), dodgeball_metalness(1.0), ray_step_size_mult(1.0);
float water_h_off(0.0), water_h_off_rel(0.0), perspective_fovy(0.0), perspective_nclip(0.0), read_mesh_zmm(0.0), indir_light_exp(1.0), cloud_height_offset(0.0);
float snow_depth(0.0), snow_random(0.0), cobj_z_bias(DEF_Z_BIAS), init_temperature(DEF_TEMPERATURE), indir_vert_offset(0.25), sm_tree_density(1.0), fog_dist_scale(1.0);
float CAMERA_RADIUS(DEF_CAMERA_RADIUS), C_STEP_HEIGHT(0.6), waypoint_sz_thresh(1.0), model3d_alpha_thresh(0.9), model3d_texture_anisotropy(1.0), dist_to_fire_sq(0.0);
float ocean_wave_height(DEF_OCEAN_WAVE_HEIGHT), tree_density_thresh(0.55), model_auto_tc_scale(0.0), model_triplanar_tc_scale(0.0), precip_dist_scale(1.0);
float custom_glaciate_exp(0.0), tree_type_rand_zone(0.0), jump_height(1.0), force_czmin(0.0), force_czmax(0.0), smap_thresh_scale(1.0), dlight_intensity_scale(1.0);
float model_mat_lod_thresh(5.0), clouds_per_tile(0.5), def_atmosphere(1.0), def_vegetation(1.0), ocean_depth_opacity_mult(1.0), erode_amount(1.0), ambient_scale(1.0);
float model_hemi_lighting_scale(0.5), pine_tree_radius_scale(1.0), sunlight_brightness(1.0), moonlight_brightness(1.0), sm_tree_scale(1.0), tt_fog_density(1.0);
float mouse_smooth_factor(0.0), tree_depth_scale(1.0);
float light_int_scale[NUM_LIGHTING_TYPES] = {1.0, 1.0, 1.0, 1.0, 1.0}, first_ray_weight[NUM_LIGHTING_TYPES] = {1.0, 1.0, 1.0, 1.0, 1.0};
double camera_zh(0.0);
point mesh_origin(all_zeros), camera_pos(all_zeros), cube_map_center(all_zeros);
string user_text, cobjs_out_fn, sphere_materials_fn, hmap_out_fn, skybox_cube_map_name, coll_damage_name, assimp_alpha_exclude_str;
colorRGB ambient_lighting_scale(1,1,1), mesh_color_scale(1,1,1);
colorRGBA flower_color(ALPHA0);
set<unsigned char> keys, keyset;
unsigned init_item_counts[] = {2, 2, 2, 6, 6}; // HEALTH, SHIELD, POWERUP, WEAPON, AMMO
vector<cube_t> smoke_bounds;

// camera variables
double c_radius(DEF_CRADIUS), c_theta(DEF_CTHETA), c_phi(DEF_CPHI), up_theta(DEF_UPTHETA), camera_y(DEF_CAMY);
float sun_rot(0.2), moon_rot(-0.2), sun_theta(1.2), moon_theta(0.3), light_factor, ball_velocity(15.0), cview_radius(1.0), player_speed(1.0);
vector3d up_vector(plus_y), cview_dir(all_zeros);
point camera_origin(all_zeros), surface_pos(all_zeros), cpos2;
char player_name[MAX_CHARS] = "Player";
bool vert_opt_flags[3] = {0}; // {enable, full_opt, verbose}


extern bool clear_landscape_vbo, use_dense_voxels, tree_4th_branches, model_calc_tan_vect, water_is_lava, use_grass_tess, def_tex_compress, ship_cube_map_reflection;
extern bool flashlight_on, player_wait_respawn, camera_in_building, player_in_tunnel;
extern int camera_flight, DISABLE_WATER, DISABLE_SCENERY, camera_invincible, onscreen_display, mesh_freq_filter, show_waypoints, last_inventory_frame;
extern int tree_coll_level, GLACIATE, UNLIMITED_WEAPONS, destroy_thresh, MAX_RUN_DIST, mesh_gen_mode, mesh_gen_shape, map_drag_x, map_drag_y, player_in_water;
extern unsigned NPTS, NRAYS, LOCAL_RAYS, GLOBAL_RAYS, DYNAMIC_RAYS, NUM_THREADS, MAX_RAY_BOUNCES, grass_density, max_unique_trees, shadow_map_sz;
extern unsigned scene_smap_vbo_invalid, spheres_mode, max_cube_map_tex_sz, DL_GRID_BS;
extern float fticks, team_damage, self_damage, player_damage, smiley_damage, smiley_speed, tree_deadness, tree_dead_prob, lm_dz_adj, nleaves_scale, flower_density, universe_ambient_scale;
extern float mesh_scale, tree_scale, mesh_height_scale, smiley_acc, hmv_scale, last_temp, grass_length, grass_width, branch_radius_scale, tree_height_scale, planet_update_rate;
extern float MESH_START_MAG, MESH_START_FREQ, MESH_MAG_MULT, MESH_FREQ_MULT, def_tex_aniso;
extern double map_x, map_y, tfticks;
extern point hmv_pos, camera_last_pos;
extern colorRGBA sunlight_color;
extern int coll_id[];
extern float tree_lod_scales[4];
extern string read_hmap_modmap_fn, write_hmap_modmap_fn, read_voxel_brush_fn, write_voxel_brush_fn, font_texture_atlas_fn;
extern vector<bbox> team_starts;
extern player_state *sstates;
extern pt_line_drawer obj_pld;
extern tree_cont_t t_trees;
extern dpart_params_t dp_params;
extern hmap_params_t hmap_params;
extern reflect_plane_selector reflect_planes;
extern reflective_cobjs_t reflective_cobjs;

// init and cleanup functions exported from other systems that are called at the beginning and end of main()
void init_keyset();
int load_config(string const &config_file);
void init_lights();

bool export_modmap(string const &filename);
void reset_planet_defaults();
void invalidate_cached_stars();
void clear_default_vao();

void create_sin_table();

void clear_sm_tree_vbos();
void clear_scenery_vbos();
void clear_asteroid_contexts();
void clear_quad_ix_buffer_context();
void clear_vbo_ring_buffer();
void free_cloud_context();
void free_universe_context();
void free_animal_context();

void setup_linear_fog(colorRGBA const &color, float fog_end);
void write_map_mode_heightmap_image();
void apply_grass_scale();
void take_screenshot_texture();
void teleport_to_map_location();
void building_gameplay_action_key(int mode, bool mouse_wheel);
float get_player_building_speed_mult();
void toggle_city_spectate_mode();

float get_tt_building_sound_gain();


// all OpenGL error handling goes through these functions
bool get_gl_error(unsigned loc_id, const char* stmt, const char* fname) {

	bool had_error(0);

	while (1) {
		int const error(glGetError());
		if (!error) break;
		const GLubyte *const error_str(gluErrorString(error));
		cerr << "GL Error " << error;
		if (fname) {cerr << " in statement " << stmt << " file " << fname << " line " << loc_id << ": ";} // from source file
		else {cerr << " at location id " << loc_id << ": ";} // from check_gl_error()
		if (error_str) {cerr << error_str << "." << endl;} else {cerr << "<NULL>." << endl;}
		had_error = 1;
	}
	if (fname) {assert(!had_error);} // caller ignores return code in this case, add the assert here
	return had_error;
}
bool check_gl_error(unsigned loc_id) {

	bool had_error(0);
#ifdef _DEBUG
	had_error = get_gl_error(loc_id);
	assert(!had_error); // currently fatal
#endif
	return had_error;
}


void display_window_resized() {invalidate_cached_stars();}
void post_window_redisplay () {glutPostRedisplay();} // Schedule a new display event


void clear_context() { // free all textures, shaders, VBOs, etc.; used on context switch and at shutdown

	reset_textures();
	free_universe_context();
	free_model_context();
	free_voxel_context();
	free_sphere_vbos();
	clear_shaders();
	reset_snow_vbos();
	update_grass_vbos();
	update_tiled_terrain_grass_vbos();
	clear_tree_context();
	clear_sm_tree_vbos();
	clear_scenery_vbos();
	reset_tiled_terrain_state();
	free_cobj_draw_group_vbos();
	clear_univ_obj_contexts();
	clear_asteroid_contexts();
	clear_quad_ix_buffer_context();
	clear_vbo_ring_buffer();
	clear_default_vao();
	free_cloud_context();
	free_animal_context();
	reflective_cobjs.free_textures();
	clear_landscape_vbo_now();
	clear_building_vbos();
	free_city_context();
}


void quit_3dworld() { // called once at the end for proper cleanup

	cout << "quitting" << endl;
	kill_current_raytrace_threads();
	end_building_rt_job();
	clear_context();
	exit_openal();

	if (!universe_only) {
		free_models();
		free_scenery_cobjs();
		delete_matrices();
	}
	glutExit();
	//_CrtDumpMemoryLeaks();
	exit(0); // quit
}


int get_swap_interval() {return (vsync_enabled ? 1 : 0);}
#ifdef _WIN32
void set_vsync() {wglSwapIntervalEXT(get_swap_interval());}
#else // linux
void set_vsync() {
  //has_extension("GLX_EXT_swap_control")
  if (glXSwapIntervalEXT) {glXSwapIntervalEXT(glXGetCurrentDisplay(), glXGetCurrentDrawable(), get_swap_interval());}
  else if (glXSwapIntervalSGI) {glXSwapIntervalSGI(get_swap_interval());}
}
#endif

void init_window() { // register all glut callbacks

	set_vsync();
	glutSetCursor(GLUT_CURSOR_CROSSHAIR);
	glutDisplayFunc(display);
 	glutReshapeFunc(resize);
	//glutCloseFunc(quit_3dworld); // can't do this because we don't want to quit when destroying the context
	// init keyboard and mouse callbacks
	glutIgnoreKeyRepeat(1);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutPassiveMotionFunc(mousePassiveMotion);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard2);
	glutKeyboardUpFunc(keyboard_up);
	glutSpecialUpFunc(keyboard2_up);
	init_keyset();

 	// Initialize GL
	fgMatrixMode(FG_PROJECTION);
 	fgLoadIdentity();
	glClearColor(0.0, 0.0, 0.0, 0.0);
}

void toggle_fullscreen() {

	static int xsz(window_width), ysz(window_height);

	if (!fullscreen) { // make fullscreen
		if (window_width  > 0) {xsz = window_width;}
		if (window_height > 0) {ysz = window_height;}
		glutFullScreen();
		glutSetCursor(GLUT_CURSOR_NONE);
	}
	else { // make windowed
		glutReshapeWindow(xsz, ysz);
		resize(xsz, ysz);
		glutSetCursor(GLUT_CURSOR_CROSSHAIR);
	}
	screen_reset = 1;
	fullscreen  ^= 1;
}


void enable_blend () {glEnable (GL_BLEND);}
void disable_blend() {glDisable(GL_BLEND);}

void set_std_blend_mode     () {glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);}
void set_additive_blend_mode() {glBlendFunc(GL_SRC_ALPHA, GL_ONE);}

void reset_fog() {
	setup_linear_fog(GRAY, ((world_mode == WMODE_INF_TERRAIN) ? get_inf_terrain_fog_dist() : 2.5*Z_SCENE_SIZE*fog_dist_scale));
}

void set_std_depth_func        () {glDepthFunc(GL_LESS  );}
void set_std_depth_func_with_eq() {glDepthFunc(GL_LEQUAL);}

void set_gl_params() {
	reset_fog();
	set_std_depth_func();
	set_std_blend_mode();
	glEnable(GL_DEPTH_TEST);
}


void reset_camera_pos() {

	if (world_mode == WMODE_UNIVERSE) return;
	camera_origin = mesh_origin;
	up_theta      = DEF_UPTHETA;
	c_radius      = DEF_CRADIUS;
	c_theta       = DEF_CTHETA;
	c_phi         = DEF_CPHI;
	camera_y      = DEF_CAMY;
	if (world_mode == WMODE_UNIVERSE) camera_origin.z = zcenter;
	if (enable_mouse_look) {m_button = GLUT_LEFT_BUTTON;}
}


void set_perspective_near_far(float near_clip, float far_clip, float aspect_ratio) {

	if (window_width == 0) return; // window not setup yet, skip (maybe got here during mouse/keyboard even before display() was called or config was read)
	assert(window_width > 0 && window_height > 0);
	perspective_nclip = near_clip;
	fgMatrixMode(FG_PROJECTION);
	fgLoadIdentity();
	fgPerspective(perspective_fovy, ((aspect_ratio == 0.0) ? ((double)window_width)/window_height : aspect_ratio), perspective_nclip, far_clip);
	fgMatrixMode(FG_MODELVIEW);
	fgLoadIdentity();
}


void set_perspective(float fovy, float nc_scale) {

	perspective_fovy = fovy;
	set_perspective_near_far(nc_scale*NEAR_CLIP, FAR_CLIP);
}


void check_zoom() {

	float fovy(PERSP_ANGLE);
	if      (do_zoom == 1) {do_zoom = 2; fovy /= ZOOM_FACTOR;}
	else if (do_zoom == 2) {do_zoom = 1;}
	set_perspective(fovy, 1.0);
}


void check_xy_offsets() {

	if (camera_view) return;
	int const mrd(((world_mode == WMODE_INF_TERRAIN) ? 4 : 1)*MAX_RUN_DIST); // increase distance in TT mode to reduce shadow map updates
	vector3d delta;
	while (xoff >=  mrd) {xoff -= mrd; delta.x -= mrd*DX_VAL;}
	while (xoff <= -mrd) {xoff += mrd; delta.x += mrd*DX_VAL;}
	while (yoff >=  mrd) {yoff -= mrd; delta.y -= mrd*DY_VAL;}
	while (yoff <= -mrd) {yoff += mrd; delta.y += mrd*DY_VAL;}
	surface_pos     += delta;
	camera_last_pos += delta;
}


float get_player_speed_mult() {return pow(3.0f, min(((world_mode == WMODE_INF_TERRAIN) ? 3 : 2), do_run));}

float calc_speed() {
	float speed(get_player_speed_mult());
	if (camera_in_air) {speed *= CAMERA_AIR_CONT;} // what about smileys?
	return speed;
}


vector3d calc_camera_direction() {return -rtp_to_xyz(1.0, c_theta, c_phi);}

void update_cpos() {

	if (world_mode != WMODE_UNIVERSE) {
		if (fabs(c_phi) < 0.001 || fabs(c_phi - PI) < 0.001 || fabs(c_phi - TWO_PI) < 0.001) c_phi += 0.01;
		if (!spectate) {cview_dir = calc_camera_direction();} // spherical coordinates
		cview_radius = c_radius;
	}
	camera_pos = camera_origin;
	if (world_mode != WMODE_UNIVERSE) {camera_pos -= vector3d_d(cview_dir)*double(cview_radius);}
}

void calc_theta_phi_from_cview_dir() {
	c_phi   = acosf(-cview_dir.z);
	c_theta = atan2(-cview_dir.y, -cview_dir.x);
}
void set_camera_pos_dir(point const &pos, vector3d const &dir) {
	//camera_pos  = pos // no, can't change this mid-frame
	surface_pos = pos;
	cview_dir   = dir.get_norm();
	calc_theta_phi_from_cview_dir();
}


void move_camera_pos_xy(vector3d const &v, float dist) {
	// called when the player moves in both ground and TT modes; checks for XY collisions that block movement
	// normal ground movement - should speed depend on orientation or not?
	if (world_mode == WMODE_INF_TERRAIN) {dist *= get_player_building_speed_mult();}
	static float prev_camera_zval(surface_pos.z); // required for walking on bridges to determine if camera is on or below the bridge
	point const prev(surface_pos);
	float const xy_scale(dist*(v.mag()/v.xy_mag()));
	surface_pos.x += xy_scale*v.x;
	surface_pos.y += xy_scale*v.y;
	if (world_mode == WMODE_INF_TERRAIN) {check_legal_movement_using_model_coll(prev, surface_pos, CAMERA_RADIUS);} // collision with models
	bool const include_cars = 1;
	proc_city_sphere_coll(surface_pos, prev, CAMERA_RADIUS, prev_camera_zval, include_cars, nullptr, 1); // use prev pos for building collisions; check_interior=1
	prev_camera_zval = surface_pos.z;
}


void move_camera_pos(vector3d const &v, float dist) { // remember that dist is negative
	if (dist == 0.0) return;
	if (!camera_surf_collide || camera_flight) {surface_pos += v*dist;}
	else {move_camera_pos_xy(v, dist);}
}

float get_player_move_dist() {return fticks*speed_mult*player_speed*GROUND_SPEED*calc_speed();}

void advance_camera(int dir) { // player movement processing

	advanced = 1;

	if (world_mode == WMODE_UNIVERSE) { // universe
		bool const hyperspeed(do_run == 2);
		if (player_ship_inited()) {player_ship().thrust(dir, speed_mult, hyperspeed);}
		return;
	}
	if (camera_mode != 1 || (map_mode && world_mode != WMODE_INF_TERRAIN)) return;
	if (world_mode == WMODE_INF_TERRAIN && player_wait_respawn) return; // can't move during respawn
	vector3d v;
	float dist(get_player_move_dist());
	
	if (game_mode && sstates != NULL) {
		if (sstates[CAMERA_ID].freeze_time > 0) return; // can't move
		dist *= sstates[CAMERA_ID].get_rspeed_scale();
	}
	// Note: the player can hold down both forward and sidestep keys at the same time and move diagonally a bit faster than they can move using forward alone;
	// while I'm sure it's possible to fix this by tracking which keys are held down each frame, I've gotten used to it, and it doesn't affect gameplay much
	switch (dir) {
	case MOVE_BACK: // backward
		dist  = -dist;
		dist *= BACKWARD_SPEED; // slower backwards
	case MOVE_FRONT: // forward
		move_camera_pos(cview_dir, dist);
		break;
	case MOVE_RIGHT:
		dist = -dist;
	case MOVE_LEFT:
		dist *= SIDESTEP_SPEED;
		cross_product(up_vector, cview_dir, v);
		move_camera_pos(v, dist);
		break;
	default: assert(0);
	}
}


void change_terrain_zoom(float val) {

	if (!(display_mode & 0x01)) { // mesh not enabled - only scale trees
		tree_scale /= val;
		regen_trees(0);
		build_cobj_tree();
		clear_tiled_terrain();
	}
	else {
		last_temp         = -100.0; // force update
		camera_change     = 1;
		mesh_scale_change = 1;
		update_mesh(val, 1);
		clear_tiled_terrain();
		calc_watershed();
	}
	scene_smap_vbo_invalid = 2; // full rebuild of shadowers
}


void change_world_mode() { // switch terrain mode: 0 = normal/ground, 1 = universe, 2 = tiled terrain

	if (map_mode || universe_only || tiled_terrain_only || (disable_universe && disable_inf_terrain)) return;
	static int xoff_(0), yoff_(0), xoff2_(0), yoff2_(0);
	static point camera_pos_(all_zeros);
	last_temp = -100.0; // force update

	if (world_mode == WMODE_UNIVERSE) { // restore saved parameters and recalculate sun and moon pos
		xoff = xoff_; yoff = yoff_; xoff2 = xoff2_; yoff2 = yoff2_;
		camera_pos = camera_pos_;
		//up_vector  = get_player_up();  // doesn't work - need to set up_theta
		//cview_dir  = get_player_dir(); // doesn't work - need to set c_theta, c_phi
		update_sun_and_moon();
	}
	else if (world_mode == WMODE_INF_TERRAIN) {
		c_radius = 2.5;
	}
	do {
		world_mode = (world_mode+1)%NUM_WMODE;
	} while ((disable_universe && world_mode == WMODE_UNIVERSE) || (disable_inf_terrain && world_mode == WMODE_INF_TERRAIN));
	
	if (world_mode == WMODE_UNIVERSE) { // save camera position parameters
		xoff_ = xoff; yoff_ = yoff; xoff2_ = xoff2; yoff2_ = yoff2;
		camera_pos_ = camera_pos; // fix_player_upv()?
	}
	else if (combined_gu) {
		setup_current_system();
	}
	if (!map_mode) {reset_offsets();} // ???
	init_x        = 1;
	camera_change = 1;
	if (world_mode == WMODE_GROUND || world_mode == WMODE_INF_TERRAIN) {apply_grass_scale();}
	reset_fog();
	clear_tiled_terrain();
	update_grass_vbos();
	clear_vbo_ring_buffer();
	obj_pld.free_mem();
	invalidate_cached_stars();
	clear_dynamic_lights();
	glDrawBuffer(GL_BACK);
	post_window_redisplay();
	if (world_mode == WMODE_GROUND && combined_gu) {regen_trees(0);}
	change_wmode_frame = frame_counter;
}


void update_sound_loops() {

	bool const universe(world_mode == WMODE_UNIVERSE);
	float const fire_gain(0.1/dist_to_fire_sq), rain_wind_volume(get_tt_building_sound_gain());
	if (player_in_tunnel) {set_sound_loop_state(SOUND_LOOP_RAIN, 1, 0.5);} // light rain sound in tunnel
	else {set_sound_loop_state(SOUND_LOOP_RAIN, (!universe && rain_wind_volume > 0.0 && is_rain_enabled()),    rain_wind_volume);}
	set_sound_loop_state(SOUND_LOOP_WIND,       (!universe && rain_wind_volume > 0.0 && wind.mag() >= 1.0),    rain_wind_volume);
	set_sound_loop_state(SOUND_LOOP_FIRE,       (!universe && dist_to_fire_sq > 0.0 && dist_to_fire_sq < 2.0), fire_gain);
	set_sound_loop_state(SOUND_LOOP_UNDERWATER, (!universe && (underwater || player_in_water == 2) && frame_counter > change_wmode_frame+1));
	dist_to_fire_sq = 0.0;
	proc_delayed_and_placed_sounds();
}


void switch_weapon(bool prev, bool mouse_wheel) {
	if (world_mode == WMODE_UNIVERSE) {player_ship().switch_weapon(prev);} else {switch_player_weapon((prev ? -1 : 1), mouse_wheel);}
}

// *** Begin glut callback functions ***

// This function is called whenever the window is resized. 
// Parameters are the new dimentions of the window
void resize(int x, int y) {

	if (init_resize) {init_resize = 0;}
	else {add_uevent_resize(x, y);}
	y = y & (~1); // make sure y is even (required for video encoding)
	x = x & (~1); // make sure x is even (required for video encoding)
 	glViewport(0, 0, x, y);
 	window_width  = x;
 	window_height = y;
	set_perspective(PERSP_ANGLE, 1.0);
	set_gl_params();
	post_window_redisplay();
	display_window_resized();
}

bool is_shift_key_pressed() {return ((glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0);}
bool is_ctrl_key_pressed () {return ((glutGetModifiers() & GLUT_ACTIVE_CTRL ) != 0);} // used for universe mode query and building gample crouch
bool is_alt_key_pressed  () {return ((glutGetModifiers() & GLUT_ACTIVE_ALT  ) != 0);}

void update_ctrl_key_pressed() {ctrl_key_pressed = is_ctrl_key_pressed();}

void update_key_ignore_ctrl(unsigned char &key) {
	if (is_ctrl_key_pressed()) {key ^= 96;} // control key modifier changes key code
}

struct player_height_mgr_t {
	double cur_height = 0.0;
	bool inited = 0;

	void next_frame() {
		if (!inited) {cur_height = camera_zh; inited = 1;} // start at full height
		float adj_val(0.0);

		if (player_wait_respawn) { // player is on the floor if waiting for respawn
			if (cur_height == 0) return; // already fully down
			adj_val = -0.5;
		}
		else if (ctrl_key_pressed) { // crouch
			if (cur_height == 0) return; // already fully down
			adj_val = -1.0;
		}
		else { // uncrouch
			if (cur_height == camera_zh) return; // already fully up
			adj_val = 1.0;
		}
		cur_height += adj_val*camera_zh*5.0*(fticks/TICKS_PER_SECOND); // 200ms for complete transition
		cur_height  = max(0.25*camera_zh, min(camera_zh, cur_height)); // clamp to valid range, 25% to 100% height
	}
};

player_height_mgr_t player_height_mgr;

double get_player_height() {return player_height_mgr.cur_height;} // control key = crouch
float  get_crouch_amt   () {return (1.0 - get_player_height()/camera_zh)/0.75;}
void force_player_height(double height) {assert(height >= 0.0); player_height_mgr.cur_height = height;} // for building attic forced crouch

// This function is called whenever the mouse is pressed or released
// button is a number 0 to 2 designating the button
// state is 1 for release 0 for press event
// x and y are the location of the mouse (in window-relative coordinates)
void mouseButton(int button, int state, int x, int y) {

	bool const fire_button(!map_mode && (button == GLUT_RIGHT_BUTTON || (enable_mouse_look && button == GLUT_LEFT_BUTTON)));
	add_uevent_mbutton(button, state, x, y);
	if (ui_intercept_mouse(button, state, x, y, 1)) return; // already handled

	if ((camera_mode == 1 || world_mode == WMODE_UNIVERSE) && fire_button) {
		b2down = !state;
		return;
	}
	m_button     = button;
	last_mouse_x = x;
	last_mouse_y = y;
	mouse_state  = state;
	if (button == GLUT_MIDDLE_BUTTON && state == 0) {user_action_key = 1;}
	
	if (state == 0) { // use mouse scroll wheel to switch weapons and zoom in/out in map mode
		if (button == 3) { // mouse wheel up
			if (map_mode) {map_zoom /= 1.2;}
			else {switch_weapon(0, 1);}
		}
		else if (button == 4) { // mouse wheel down
			if (map_mode) {map_zoom *= 1.2;}
			else {switch_weapon(1, 1);}
		}
	}
}


void clamp_and_scale_mouse_delta(int &delta, int dmax) {
	delta = min(dmax, max(-dmax, delta));
	if (abs(delta) > 1) {delta = round_fp(mouse_sensitivity*delta);} // don't round +/-1 down to 0
}

void update_mouse_pos(int dx, int dy) {

	if (animate && mouse_smooth_factor > 0.0) { // mouse smoothing
		// Note: this doesn't work well because it's only run when the mouse moves rather than every frame
		static float prev_tfticks(0.0);
		static float prev_dx(0), prev_dy(0);
		if (dx == 0 && dy == 0 && tfticks == prev_tfticks) return; // same frame, no change
		float const elapsed_time(tfticks - prev_tfticks), blend(min(1.0f*elapsed_time/mouse_smooth_factor, 1.0f));
		prev_dx = blend*dx + (1.0 - blend)*prev_dx;
		prev_dy = blend*dy + (1.0 - blend)*prev_dy;
		int const new_dx(round_fp(prev_dx)), new_dy(round_fp(prev_dy));
		if (new_dx != 0) {dx = new_dx;} else if (dx != 0) {dx = ((dx < 0) ? -1 : 1);}
		if (new_dy != 0) {dy = new_dy;} else if (dy != 0) {dy = ((dy < 0) ? -1 : 1);}
		prev_tfticks = tfticks;
	}
	int button(m_button);
	if (camera_mode == 1 && enable_mouse_look && !map_mode) {button = GLUT_LEFT_BUTTON;}

	switch (button) {
	case GLUT_LEFT_BUTTON: // h: longitude, v: latitude
		if (dx == 0 && dy == 0) break;

		if (world_mode == WMODE_UNIVERSE) {
			vector3d delta(dx, dy, 0.0); // mouse delta
			delta *= (1280.0/window_width); // ???
			if (player_ship_inited()) {player_ship().turn(delta);}
			update_cpos();
		}
		else if (map_mode) { // map mode click and drag
			if (mouse_state == 0) {map_drag_x -= dx; map_drag_y += dy;} // mouse down
		}
		else { // ground mode
			float c_phi2(c_phi - double(MOUSE_ANG_ADJ)*dy);
			if (camera_mode) {c_phi2 = max(0.01f, min((float)PI-0.01f, c_phi2));} // walking on ground

			if (dy > 0) { // change camera y direction when camera moved through poles and jump over poles (x=0,z=0) to eliminate "singularity"
				if (c_phi2 < 0.0 || (c_phi > PI && c_phi2 < PI)) {camera_y *= -1.0;}
				if (fabs(c_phi2) < MA_TOLERANCE || fabs(c_phi2 - PI) < MA_TOLERANCE) {++dy;}
			}
			else if (dy < 0) {
				if (c_phi2 > TWO_PI || (c_phi < PI && c_phi2 > PI)) {camera_y *= -1.0;}
				if (fabs(c_phi2 - TWO_PI) < MA_TOLERANCE || fabs(c_phi2 - PI) < MA_TOLERANCE) {--dy;}
			}
			c_theta -= double(MOUSE_ANG_ADJ)*dx*camera_y;
			c_theta  = fix_angle(c_theta);
			c_phi    = fix_angle(c_phi2);
			update_cpos();
		}
		break;

	case GLUT_MIDDLE_BUTTON: // translate camera
		if (!camera_view) {
			camera_origin.x += MOUSE_TRAN_ADJ*dy;
			camera_origin.y += MOUSE_TRAN_ADJ*dx;
		}
		// Note: could use the middle button to move the sun/moon, etc.
		break;

	case GLUT_RIGHT_BUTTON: // v: radius, h: up_vector
		if (map_mode) { // map scroll
			if      (dy < 0) {map_zoom /= (1.0 - 0.01*dy);}
			else if (dy > 0) {map_zoom *= (1.0 + 0.01*dy);}
			break;
		}
		if (camera_mode == 1) break; // not free look mode

		if (!camera_view) {
			up_theta += double(MOUSE_ANG_ADJ)*dx;
			up_theta  = fix_angle(up_theta);
			c_radius  = c_radius*(1.0 + double(MOUSE_R_ADJ)*dy);
			if (c_radius <= 0.05*MOUSE_R_ADJ) {c_radius = 0.05*MOUSE_R_ADJ;}
			update_cpos();
		}
		break;

	default: // is there any other mouse button? error?
		break;
	}
}

// This function is called whenever the mouse is moved with a mouse button held down.
// x and y are the location of the mouse (in window-relative coordinates)
void mouseMotion(int x, int y) {

	if (screen_reset || start_maximized) {
		last_mouse_x = x;
		last_mouse_y = y;
		screen_reset = 0;
		return;
	}
	add_uevent_mmotion(x, y);
	if (ui_intercept_mouse(0, 0, x, y, 0)) return; // already handled
	int dx(x - last_mouse_x), dy(y - last_mouse_y);
	clamp_and_scale_mouse_delta(dx, window_width /20); // limit to a reasonable delta in case the frame rate is very low
	clamp_and_scale_mouse_delta(dy, window_height/20);
	update_mouse_pos(dx, dy);
	last_mouse_x = x;
	last_mouse_y = y;

	if (enable_mouse_look) { // wrap the pointer when it goes off screen
		if (x == 0) {
			glutWarpPointer(window_width-2, y);
			last_mouse_x = window_width-2;
		}
		if (x >= window_width-1) {
			glutWarpPointer(1, y);
			last_mouse_x = 1;
		}
		if (y == 0) {
			glutWarpPointer(x, window_height-2);
			last_mouse_y = window_height-2;
		}
		if (y >= window_height-1) {
			glutWarpPointer(x, 1);
			last_mouse_y = 1;
		}
	}
	if (dx != 0 || dy != 0) {post_window_redisplay();}
}


void mousePassiveMotion(int x, int y) {
	if (enable_mouse_look) {mouseMotion(x, y);}
}

void change_tree_mode() {

	if (world_mode != WMODE_GROUND && world_mode != WMODE_INF_TERRAIN) return;
	if (num_trees == 0 && t_trees.empty()) return;
	tree_mode = (tree_mode+1)%4; // 0=none, 1=large, 2=small, 3=large+small
			
	if (world_mode == WMODE_INF_TERRAIN) {
		clear_tiled_terrain(1); // no_regen_buildings=1
	}
	else {
		//if (num_trees == 0) return; // Note: will skip scene/cobj updates on scenes that have placed trees
#if 1
		gen_scene(0, 1, 1, 0, 1); // Note: will destroy any fixed cobjs
#else
		remove_small_tree_cobjs();
		remove_tree_cobjs();
		regen_trees(0); // Note: won't regen trees if num_trees == 0
#endif
	}
}

void switch_weapon_mode() {
	if (sstates == NULL || game_mode != GAME_MODE_FPS) return;
	++sstates[CAMERA_ID].wmode;
	sstates[CAMERA_ID].verify_wmode();
	play_switch_wmode_sound();
	//last_inventory_frame = frame_counter;
}

void toggle_camera_mode() {
	camera_mode   = (camera_mode == 0);
	camera_reset  = 1;
	camera_change = 1;
	if (camera_mode == 1) {camera_invincible = 1;} // in air (else on ground)
}

void update_precip_rate_verbose(float val) {
	update_precip_rate(val);
	cout << ((val > 1.0) ? "increase" : "decrease") << " precip to " << obj_groups[coll_id[PRECIP]].max_objs << endl;
}
void show_bool_option_change(string const &name, bool new_val) {
	print_text_onscreen((name + (new_val ? " ON" : " OFF")), WHITE, 1.0, 1.0*TICKS_PER_SECOND);
}
void show_speed() {
	ostringstream oss;
	oss << "Player Speed " << get_player_speed_mult() << "x";
	print_text_onscreen(oss.str(), WHITE, 1.0, 1.0*TICKS_PER_SECOND);
}
void show_text_prompt() {
	print_text_onscreen("Text Prompt", WHITE, 1.0, 2.0*TICKS_PER_SECOND);
}

void print_texture_stats();
void print_shader_stats();
void show_tiled_terrain_debug_stats();

void show_frame_stats() {
	print_texture_stats();
	print_shader_stats();
	cout << "Draw calls for frame " << frame_counter << ": " << num_frame_draw_calls << endl;
	if (world_mode == WMODE_INF_TERRAIN) {show_tiled_terrain_debug_stats();}
}

void next_game_mode() {
	if (world_mode == WMODE_UNIVERSE)      return; // only one game mode
	if ((game_mode_disable_mask & 7) == 7) return; // all game modes disabled, leave at init game mode (error?)
	int const prev_game_mode(game_mode);
	do {game_mode = (game_mode + 1) % 3;} while (game_mode_disable_mask & (1 << game_mode)); // select the next enabled game mode
	if (game_mode != prev_game_mode) {change_game_mode();}
}


// This function is called whenever there is a keyboard input;
// key is the ASCII value of the key pressed (esc = 27, enter = 13, backspace = 8, tab = 9, del = 127);
// x and y are the location of the mouse, which generally aren't used but are part of the callback function;
// key repeat is only enabled for movement keys wasd
void keyboard_proc(unsigned char key, int x, int y) {

    switch (key) { // available: sometimes O, sometimes Z
	case 0x1B: // ESC key (27)
		quit_3dworld();
		break;
	case 'Q':
		reload_all_shaders();
		break;

	case 'A':
		enable_multisample ^= 1;
		show_bool_option_change("Multisample", enable_multisample);
		if (!enable_multisample) {glDisable(GL_MULTISAMPLE);}
		break;

	case 'X': // change selected UI menu
		next_selected_menu_ix();
		break;

	case 8: // backspace
		if (world_mode == WMODE_INF_TERRAIN) {inf_terrain_undo_hmap_mod();}
		else if (world_mode == WMODE_GROUND) {undo_voxel_brush();}
		break;
	
	case 'm': // toggle fullscreen mode
		toggle_fullscreen();
		break;

	case 'r': // run mode (always move forward)
		if (world_mode == WMODE_INF_TERRAIN && have_buildings()) {building_gameplay_action_key(3, 0);}
		else {run_forward = !run_forward;}
		break;
	case 'V': // change mouse mode
		enable_mouse_look = !enable_mouse_look;
		break;

	case 'C': // recreate mesh / add red ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_RED); // red
			break;
		}
		if (mesh_seed != 0 || read_heightmap) break;
		rand_gen_index = mesh_rgen_index = rand();
		gen_scene(1, (world_mode == WMODE_GROUND), 0, 0, 0);
		regen_lightmap();
		recreated     = 1;
		camera_change = 1;
		break;

	case 'x': // toggle auto framerate updates
		animate = !animate;
		if (animate) {reset_timing = 1;}
		break;
	case 't': // animation/physics timesteps - movement (freeze frame on objects), show star streams in universe mode
		animate2 = !animate2;
		show_bool_option_change("Timestep Update", animate2);
		if (animate2) {reset_timing = 1;}
		break;

	case 'b': // begin motion animation
		begin_motion = !begin_motion;
		show_bool_option_change("Physics Update", begin_motion);
		free_dodgeballs(1, 1);
		break;

	case 'y': // save eventlist
		save_ueventlist();
		break;

	case 'k': // change drawing model for mesh (filled polygons vs. wireframe)
		if (map_mode) {map_color = !map_color;} else {draw_model = !draw_model;}
		break;

	case 'i': // toggle autopilot / change pedestrian animation / force reflection cube map update
		if (world_mode == WMODE_UNIVERSE) {toggle_autopilot();}
		else if (world_mode == WMODE_INF_TERRAIN) {next_pedestrian_animation();}
		else {force_ref_cmap_update ^= 1;}
		break;

	case 'n': // toggle fog / reset player target
		if (world_mode == WMODE_UNIVERSE) {
			reset_player_target();
			break;
		}
		show_fog = !show_fog;
		break;

	case 'l': // enable lightning / dock fighters
		if (world_mode == WMODE_UNIVERSE) {
			toggle_dock_fighters();
			break;
		}
		show_lightning = !show_lightning;
		break;

	case 'z': // zoom / onscreen display
		if (map_mode) {map_zoom /= MAP_ZOOM;}
		else if (world_mode == WMODE_UNIVERSE) {onscreen_display = !onscreen_display;}
		else if (world_mode == WMODE_GROUND || world_mode == WMODE_INF_TERRAIN) {do_zoom = !do_zoom;} // Note: tiled mesh reflection is incompatible with zoom and is disabled
		break;

	case 'f': // print framerate and stats
		show_framerate = 1;
		timing_profiler_stats();
		show_frame_stats();
		break;
	case 'g': // pause/resume playback of eventlist
		show_bool_option_change("Frame Pause", !pause_frame);
		pause_frame = !pause_frame;
		break;
	case 'G': // toggle show framerate/universe stats / voxel add/remove (used to be z)
		stats_display_mode = (stats_display_mode + 1) % 4;
		display_framerate  = (stats_display_mode != 3);
		break;

	case 'p': // reset camera / change fire primary
		if (world_mode == WMODE_UNIVERSE) {
			change_fire_primary();
			break;
		}
		reset_camera_pos();
		break;

	case 'P': // voxel add/remove toggle
		voxel_add_remove ^= 1;
		break;

	case 'v': // reset camera and change camera mode from air to surface
		reset_camera_pos();
		if (world_mode == WMODE_INF_TERRAIN) {camera_reset = camera_change = 1; break;} // reset camera in case player died
		if (world_mode != WMODE_GROUND) break; // universe/inf terrain mode
		gamemode_rand_appear();
		toggle_camera_mode();
		break;

	case 'h': // change camera surface collision detection
		if (world_mode == WMODE_UNIVERSE) {claim_planet = 1; break;} // player claim nearby planet
		camera_surf_collide = !camera_surf_collide;
		camera_change       = 1;
		// reset last_pos so that the camera doesn't snap back to the old pos when clipping is re-enabled
		if (camera_surf_collide) {camera_last_pos = surface_pos;}
		if (world_mode == WMODE_INF_TERRAIN) {show_bool_option_change("Player Collision", camera_surf_collide);}
		break;

	case 'j': // smooth camera collision detection / hold fighters / teleport to screenshot
		if (world_mode == WMODE_UNIVERSE) {toggle_hold_fighters();}
		else if (world_mode == WMODE_INF_TERRAIN) {teleport_to_screenshot = 1;}
		else {camera_coll_smooth = !camera_coll_smooth;} // ground mode
		break;

		// camera movement
	case 'w': // advance surface camera forward
		advance_camera(MOVE_FRONT); break;
	case 's': // advance surface camera backwards
		advance_camera(MOVE_BACK);  break;
	case 'a': // step left
		advance_camera(MOVE_LEFT);  break;
	case 'd': // step right
		advance_camera(MOVE_RIGHT); break;

	case 'q': // previous weapon
		switch_weapon(1, 0);
		break;
	case 'e': // next weapon
		switch_weapon(0, 0);
		break;
	case 'W': // switch weapon mode
		switch_weapon_mode();
		break;

	case 'o': // toggle vsync
		vsync_enabled ^= 1;
		show_bool_option_change("Vsync", vsync_enabled);
		set_vsync();
		break;
	case 'u': // toggle timing profiler
		toggle_timing_profiler(); // show_bool_option_change()?
		break;
	case 'O': // text prompt (for future use)
		show_text_prompt();
		break;

	case '=': // increase temp
		temperature += TEMP_INCREMENT;
		cout << "Temperature = " << temperature << " degrees C." << endl;
		temp_change = 1;
		break;
	case '-': // decrease temp
		temperature -= TEMP_INCREMENT;
		temperature  = max(temperature, ABSOLUTE_ZERO);
		cout << "Temperature = " << temperature << " degrees C." << endl;
		temp_change = 1;
		break;

	case '\'': // increase timestep
		change_timestep(D_TIMESTEP);     break;
	case ';': // decrease timestep
		change_timestep(1.0/D_TIMESTEP); break;

	case 'T': // delete and generate tree(s) and scenery
		if (world_mode == WMODE_UNIVERSE) {
			player_ship().reset_ammo();
		}
		else {
			rand_gen_index = rand(); // Note: doesn't set mesh_rgen_index

			if (world_mode == WMODE_GROUND) {
				//gen_scenery(t_trees);
				//regen_trees(0);
				gen_scene(0, 1, 1, 0, 1);
			}
			else {
				clear_tiled_terrain();
			}
		}
		break;

	case 'R': // run mode
		++do_run;
		if (do_run > 3) {do_run = 0;}
		if (world_mode == WMODE_UNIVERSE) {change_speed_mode(do_run);}
		else {show_speed();} // onscreen printout
		//cout << "run mode = " << do_run << endl;
		break;

	case 'K': // toggle overhead map mode / ship reflections
		if (world_mode == WMODE_UNIVERSE) {ship_cube_map_reflection ^= 1;}
		else if (map_mode) {map_mode = 0;} else {map_mode = 2;}
		break;

	case 'U':
		if (/*!disable_universe &&*/ world_mode != WMODE_UNIVERSE) { // toggle universe background mode
			combined_gu = !combined_gu;
			
			if (combined_gu) { // do a fake draw pass to force the universe to be created so we can determine the closest planet/moon and setup lighting/water/temperature/vegetation/etc.
				draw_universe(1, 1, 1, 2, 1, 1); // gen_only=1
				setup_current_system();
			}
			else {
				reset_planet_defaults(); // have to do this so that regen_trees gets correct vegetation
			}
			if (world_mode == WMODE_GROUND) { // not TT
				remove_tree_cobjs();
				regen_trees(0);
				build_cobj_tree();
				gen_grass();
			}
			clear_tiled_terrain(1); // no_regen_buildings=1
			
			if (!combined_gu) {
				setup_landscape_tex_colors(ALPHA0, ALPHA0);
				disable_light(get_universe_ambient_light(1)); // disable universe ambient (not required?)
				calc_visibility(SUN_SHADOW); // reclaculate sun
				DISABLE_WATER = INIT_DISABLE_WATER;
			}
			create_landscape_texture();
		}
		break;

	case 'S': // ship stop / fog mode / toggle TT room lights
		if (world_mode == WMODE_UNIVERSE) {toggle_player_ship_stop(); break;}
		else if (world_mode == WMODE_GROUND) {use_smoke_for_fog = (use_smoke_for_fog+1) % 3;} // {normal smoke, smoke with noise, fog as smoke}
		else if (world_mode == WMODE_INF_TERRAIN) {toggle_room_light = 1;}
		break;
	case 'Z':
		if (map_mode) {map_zoom *= MAP_ZOOM; break;}
		else if (world_mode == WMODE_INF_TERRAIN) {building_action_key = 1;}
		// else avialable
		break;

	case 'E': // reset leaf colors / building gameplay use object
		if (world_mode == WMODE_INF_TERRAIN && have_buildings()) {building_gameplay_action_key(2, 0); break;} // use building object
		leaf_color_coherence = 0.5;
		tree_color_coherence = 0.2;
		leaf_base_color.R    = 0.2;
		leaf_base_color.G    = 1.0;
		register_leaf_color_change();
		break;

	case 'L': // increase terrain zoom
		if (mesh_seed != 0 || read_heightmap) break;
		if (world_mode != WMODE_UNIVERSE) {change_terrain_zoom(2.0);}
		break;
	case 'Y': // decrease terrain zoom / select closest target ship
		if (world_mode == WMODE_UNIVERSE) {auto_target_player_closest_enemy();}
		else if (!(mesh_seed != 0 || read_heightmap)) {change_terrain_zoom(0.5);}
		break;

	// object enables
	case 'c': // toggle smileys / add blue ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_BLUE); // blue
			break;
		}
		free_dodgeballs(0, 1);
		obj_groups[coll_id[SMILEY]].toggle_enable();
		if (obj_groups[coll_id[SMILEY]].enabled) {init_smileys();}
		break;
	case 'B': // precipitation / add neutral ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_NEUTRAL); // neutral
			break;
		}
		precip_mode = (precip_mode + 1) & ((world_mode == WMODE_INF_TERRAIN) ? 1 : 3); // 4 modes: 0=none, 1=new, 2=old, 3=new+old
		is_cloudy   = (precip_mode > 0);
		if (!(precip_mode & 1)) {obj_groups[coll_id[PRECIP]].toggle_enable();}
		//if (obj_groups[coll_id[PRECIP]].is_enabled()) {seed_water_on_mesh(10.0);} // instantly seed the mesh with water
		break;

	case 'N': // decrease precipitation rate by 1.5X
		update_precip_rate_verbose(1.0/1.5);
		break;
	case 'M': // increase precipitation rate by 1.5X
		update_precip_rate_verbose(1.5);
		break;

	case 'H': // save mesh state/modmap/voxel brushes/cobj file/materials/heightmap
		if (world_mode == WMODE_UNIVERSE) {export_modmap("output.modmap");}
		else if (map_mode) {write_map_mode_heightmap_image();}
		else if (world_mode == WMODE_GROUND) {
			if (voxel_editing) {write_voxel_brushes();}
			else if (spheres_mode) {
				if (show_scores) {write_sphere_materials_file(sphere_materials_fn);} // okay if fails
				else {write_def_coll_objects_file();}
			}
			else {save_state(state_file);}
		}
		else if (world_mode == WMODE_INF_TERRAIN) {write_default_hmap_modmap();}
		break;
	case 'J': // load mesh state / take TT screenshot
		if (world_mode == WMODE_GROUND) {load_state(state_file);}
		else if (world_mode == WMODE_INF_TERRAIN) {take_screenshot_texture();}
		break;
	case 'I': // write mesh points / toggle building interiors
		if (world_mode == WMODE_GROUND) {write_mesh(mesh_file);}
		else if (world_mode == WMODE_INF_TERRAIN) {draw_building_interiors ^= 1;}
		break;

	// screenshots
	case 'D': // .bmp
		screenshot(window_width, window_height, "./", 1);
		break;
	case 'F': // .jpg
		screenshot(window_width, window_height, "./", 0);
		break;

	// rotate sun/moon
	case '[':
		sun_rot  += LIGHT_ROT_AMT;
		update_sun_shadows();
		break;
	case ']':
		sun_rot  -= LIGHT_ROT_AMT;
		update_sun_shadows();
		break;
	case '{':
		moon_rot += LIGHT_ROT_AMT;
		calc_visibility(MOON_SHADOW);
		break;
	case '}':
		moon_rot -= LIGHT_ROT_AMT;
		calc_visibility(MOON_SHADOW);
		break;

	case ' ': // fire/jump/respawn key
		if (world_mode == WMODE_GROUND && camera_mode == 1 && camera_surf_collide && enable_mouse_look) { // jump
			if (!spectate && sstates != nullptr) {sstates[CAMERA_ID].jump(get_camera_pos());}
		}
		else if (world_mode == WMODE_GROUND && game_mode && camera_mode == 0 && !spectate && sstates != nullptr && sstates[CAMERA_ID].deaths > 0 && !(sstates[CAMERA_ID].wmode&1)) { // respawn
			gamemode_rand_appear();
			toggle_camera_mode();
			sstates[CAMERA_ID].jump_time = 0.25*TICKS_PER_SECOND; // suppress extra jump if space is held down too long
		}
		else if (world_mode == WMODE_INF_TERRAIN && map_mode) {teleport_to_map_location();}
		else {fire_weapon(); flashlight_on = 0;} // fire; spacebar doesn't toggle flashlight because flashlight_on is cleared at the beginning of display()
		break;
	case '<': // decrease weapon velocity
		ball_velocity = max(0.0, ball_velocity-5.0);
		break;
	case '>': // decrease weapon velocity
		ball_velocity += 5.0;
		break;

	case '	': // tab
		show_scores = !show_scores;
		break;

	case '1': // toggle mesh draw / universe star/planet/moon distance culling
		display_mode ^= 0x01;   break;
	case '2': // toggle grass/snow draw / universe colonization coloring
		display_mode ^= 0x02;   break;
	case '3': // toggle water/ice
		display_mode ^= 0x04;   break;
	case '4': // toggle occlusion culling / tiled terrain/voxel/mesh detail normal maps
		display_mode ^= 0x08;   show_bool_option_change("Occlusion Culling", (display_mode & 0x08)); break;
	case '5': // walk on snow/ship shadows/reflections/debugging
		display_mode ^= 0x10;   show_bool_option_change((camera_in_building ? "Indirect Lighting" : "Debug Mode"), (display_mode & 0x10)); break;
	case '6': // toggle water reflections, bump maps, bloom, and map view lighting/shadows
		display_mode ^= 0x20;   break;
	case '7': // toggle snow accumulation, clouds, and universe mode multithreading
		display_mode ^= 0x40;   break;
	case '8': // toggle water caustics/smoke accumulation and DOF
		display_mode ^= 0x80;   break;
	case '9': // toggle leaf wind, ocean waves, footsteps, snow footprints, asteroid belt fog, flashlight indirect, and dynamic particle drawing
		display_mode ^= 0x0100; break;
	case '0': // toggle universe stencil shadows / toggle spraypaint mode / toggle particles / toggle TT tree leaf shadows
		if (world_mode == WMODE_UNIVERSE) {univ_stencil_shadows ^= 1;}
		else if (world_mode == WMODE_GROUND) {
			if (begin_motion || show_scores) {toggle_sphere_mode();} else {toggle_spraypaint_mode();}
		}
		else {display_mode ^= 0x0200;}
		break;

	case '\\': // enable dynamic particles (to test dynamic lighting, dynamic shadows, and collision detection) / helicopter shadow updates
		if (world_mode == WMODE_INF_TERRAIN) {enable_hcopter_shadows ^= 1; break;}
		display_mode ^= 0x0200;
		d_part_sys.clear();
		break;
	}
	post_window_redisplay();
}


void print_wind() {cout << "wind: " << wind.str() << endl;}
double get_map_shift_val() {return map_zoom*double(MAP_SHIFT)*(is_shift_key_pressed() ? 8 : 1);}


// handles user key remapping and disabling of keys in gameplay mode
class keyboard_remap_t {

	map<int, int> key_map;
	set<int> key_null, enabled_keys, disabled_keys;

	static int get_numeric_char(char c) {
		if (c >= '0' && c <= '9') {return (c - '0');}
		if (c >= 'a' && c <= 'f') {return (c - 'a' + 10);}
		if (c >= 'A' && c <= 'F') {return (c - 'A' + 10);}
		cerr << "Error extracting hex value from character '" << c << "'" << endl;
		return -1;
	}
	static int extract_char_or_hex_number(string const &str) {
		if (str.size() == 1) { // single character
			return int(str[0]);
		}
		if (str.size() >= 3 && str[0] == '0' && (str[1] == 'x' || str[1] == 'X')) { // hex number
			if (str.size() == 3) { // 1-digit
				return get_numeric_char(str[2]);
			}
			else if (str.size() == 4) { // 2-digit
				int const hi(get_numeric_char(str[2])), lo(get_numeric_char(str[3]));
				return ((hi < 0 || lo < 0) ? -1 : (16*hi + lo));
			}
		}
		cerr << "Error extracting char or hex number from string '" << str << "'" << endl;
		return -1;
	}
	static void add_keys_to_set(string const &keys_to_add, set<int> &key_set) { // char vs. int?
		for (string::const_iterator i = keys_to_add.begin(); i != keys_to_add.end(); ++i) {key_set.insert(*i);}
	}

public:
	void add_key_null(int key) {
		assert(key >= 0);
		key_null.insert(key);
		key_map.erase(key); // just in case
	}
	void add_key_remap(int key_from, int key_to) {
		assert(key_from >= 0 && key_to >= 0);
		// okay to map a key to itself, remap a key that was already remapped, or remap multiple keys to the same value
		key_map[key_from] = key_to;
	}
	bool parse_remap_command(string const &sfrom, string const &sto) {
		assert(!sfrom.empty() && !sto.empty());
		int const kfrom(extract_char_or_hex_number(sfrom));
		if (kfrom < 0) return 0;

		if (sto == "null" || sto == "NULL") {
			add_key_null(kfrom);
		}
		else {
			int const kto(extract_char_or_hex_number(sto));
			if (kto < 0) return 0;
			add_key_remap(kfrom, kto);
		}
		return 1;
	}
	template<typename key_t> bool remap_key(key_t &key, bool special, bool up) const { // unsigned char and int (up is ignored)
		if (special) return 1; // special keys can't be remapped yet
		if (key_null.find(key) != key_null.end()) return 0; // null key, ignore
		auto it(key_map.find(key));
		if (it != key_map.end()) {key = key_t(it->second);}
		return is_key_enabled(key);
	}

	// key enabling/disabling code (doesn't apply to special or modifier keys)
	bool is_key_enabled(int key) const {
		if (key == 0x1B) return 1; // we always enable the escape/quit key
		if (disabled_keys.find(key) != disabled_keys.end()) return 0; // key explicitly disabled
		if (enabled_keys.empty()) return 1; // no enabled keys => all keys are enabled
		return (enabled_keys.find(key) != enabled_keys.end()); // check enabled set
	}
	void enable_only_keys(string const &ekeys) { // empty keys = all
		enabled_keys.clear();
		add_keys_to_set(ekeys, enabled_keys);
	}
	void set_disabled_keys(string const &dkeys) {
		disabled_keys.clear();
		add_keys_to_set(dkeys, disabled_keys);
	}
};

keyboard_remap_t kbd_remap;


void keyboard2(int key, int x, int y) { // handling of special keys

	if (ui_intercept_keyboard(key, 1))   return; // already handled
	if (!kbd_remap.remap_key(key, 1, 0)) return;
	add_uevent_keyboard_special(key, x, y);
#ifdef _WIN32
	update_ctrl_key_pressed();
#else // linux
	if (key == 114) {ctrl_key_pressed = 1;} // key code 114 appears to be the CTRL key
#endif

	switch (key) { // unused: F9
	case GLUT_KEY_UP:
	case GLUT_KEY_DOWN:
		if (map_mode) {map_y += ((key == GLUT_KEY_UP) ? 1 : -1)*(get_map_shift_val() + 0);}
		else {
			wind.y += ((key == GLUT_KEY_UP) ? 1 : -1)*WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_LEFT:
	case GLUT_KEY_RIGHT:
		if (map_mode) {map_x += ((key == GLUT_KEY_RIGHT) ? 1 : -1)*(get_map_shift_val() + 0);}
		else {
			wind.x += ((key == GLUT_KEY_RIGHT) ? 1 : -1)*WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_F1: // switch terrain mode: 0 = normal/ground, 1 = universe, 2 = infinite terrain
		change_world_mode();
		break;

	case GLUT_KEY_F2: // switch game mode
		next_game_mode();
		break;

	case GLUT_KEY_F3:
		enable_postproc_recolor ^= 1;
		break;

	case GLUT_KEY_F4: // switch weapon mode
		if (world_mode != WMODE_UNIVERSE) {switch_weapon_mode();}
		break;
	case GLUT_KEY_F5: // toggle large/small trees
		change_tree_mode();
		break;

	case GLUT_KEY_F6: // enable/disable gameplay mode keys
		{
			static bool gameplay_key_mode(0);
			gameplay_key_mode ^= 1;
			// empty string enables all keys when not in gameplay_key_mode
			kbd_remap.enable_only_keys(gameplay_key_mode ? "asdwqe " : "");
			print_text_onscreen((gameplay_key_mode ? "Disabling Non-Gameplay Keys" : "Enabling All Keys"), PURPLE, 1.0, TICKS_PER_SECOND, 10);
		}
		break;

	case GLUT_KEY_F7: // toggle auto time advance
		obj_groups[coll_id[PRECIP]].app_rate = 40;
		auto_time_adv = (auto_time_adv+1)%5;
		cout << "Auto time advance = " << auto_time_adv << endl;
		break;

	case GLUT_KEY_F8: // toggle spectator gameplay mode
		if (world_mode == WMODE_INF_TERRAIN && have_buildings()) {toggle_city_spectate_mode(); break;}
		if (!spectate && (num_smileys == 0 || !obj_groups[coll_id[SMILEY]].enabled)) break;
		if (spectate) {camera_reset = camera_change = 1;}
		spectate = !spectate;
		show_bool_option_change("Spectate", spectate);
		break;

	case GLUT_KEY_F9: // unused
		break;

	case GLUT_KEY_F10: // switch cloud model / toggle smoke_dlights
		cloud_model = !cloud_model;
		smoke_dlights ^= 1;
		show_bool_option_change("Smoke Dynamic Lights", smoke_dlights);
		break;
	case GLUT_KEY_F11: // temporary toggle of core context mode
		use_core_context ^= 1;
		show_bool_option_change("Core Context", use_core_context);
		break;
	case GLUT_KEY_F12: // toggle volumetric lighting
		volume_lighting ^= 1;
		show_bool_option_change("Volume Lighting", volume_lighting);
		break;
	}
	post_window_redisplay();
}


void init_keyset() {

	string keyvals = "wsad "; // movement; these keys repeat each frame
	for (unsigned i = 0; i < keyvals.size(); ++i) {keyset.insert(keyvals[i]);}
}


unsigned char get_key_other_case(unsigned char key) {

	unsigned char key2(key); // check if shift key was pressed while holding down a key
	if (key >= 'a' && key <= 'z') key2 = key + ('A' - 'a');
	if (key >= 'A' && key <= 'Z') key2 = key + ('a' - 'A');
	return key2;
}


void keyboard_up(unsigned char key, int x, int y) {

	update_key_ignore_ctrl(key);
	if (kbd_text_mode || key == 13) return; // ignore text mode and enter key
	if (!kbd_remap.remap_key(key, 0, 1)) return;
	if (keyset.find(key) != keyset.end()) {add_uevent_keyboard_up(key, x, y);}
	keyset_it it(keys.find(key));

	if (it == keys.end()) {
		unsigned char key2(get_key_other_case(key));
		if (key2 != key) {it = keys.find(key2);}

		if (it == keys.end()) {
			if (key == 9) return; // alt-tab
			cout << "Warning: Keyboard up event for key " << key << " (" << int(key) << ") with no corresponding keyboard down event." << endl;
			//assert(0); // too strong?
			return;
		}
	}
	keys.erase(it);
}


void keyboard2_up(int key, int x, int y) {

	if (!kbd_remap.remap_key(key, 1, 1)) return;
#ifdef _WIN32
	update_ctrl_key_pressed();
#else // linux
	if (key == 114) {ctrl_key_pressed = 0;} // key code 114 appears to be the CTRL key
#endif
	// nothing else to do here
}


void exec_text(string const &text) {

	if (world_mode == WMODE_UNIVERSE) { // handled in universe_control.cpp
		exec_universe_text(text);
		return;
	}
	cout << "Text: " << text << endl;
	print_text_onscreen(text, WHITE, 1.0, TICKS_PER_SECOND, 999);
}


void keyboard(unsigned char key, int x, int y) {

	update_key_ignore_ctrl(key);
	if (ui_intercept_keyboard(key, 0)) return; // already handled (should this go into keyboard_proc()?)

	if (key == 13) { // enter key - toggle text mode
		if (kbd_text_mode) {exec_text(user_text);}
		user_text.clear();
		kbd_text_mode = !kbd_text_mode;
		if (!kbd_text_mode) return;
	}
	if (kbd_text_mode) {
		if (key == 8) { // delete key
			if (!user_text.empty()) {user_text.erase(user_text.begin()+user_text.size()-1);} // pop_back() for string
			return;
		}
		if (key != 13) {user_text.push_back(key);} // not the enter key from above (fallthrough case)
		print_text_onscreen((string("Enter Text: ") + user_text), WHITE, 1.0, 10*TICKS_PER_SECOND, 999);
		return;
	}
	if (!kbd_remap.remap_key(key, 0, 0)) return;
	add_uevent_keyboard(key, x, y);

	if (keys.find(key) != keys.end()) { // can happen with control clicks
		cout << "Warning: Keyboard event for key " << key << " (" << int(key) << ") which has alredy been pressed." << endl;
		return;
	}
	keys.insert(key);
	if (keyset.find(key) == keyset.end()) {keyboard_proc(key, x, y);}
}


void proc_kbd_events() {

	for (keyset_it it = keys.begin(); it != keys.end(); ++it) {
		if (keyset.find(*it) != keyset.end()) {keyboard_proc(*it, 0, 0);} // x and y = ?
	}
	player_height_mgr.next_frame(); // updated based on control key and run once at the beginning of each frame
	if (animate && mouse_smooth_factor > 0.0) {update_mouse_pos(0, 0);}
}


void alloc_if_req(char *&fn, const char *def_fn=nullptr) {
	if (fn == nullptr) {
		fn = new char[MAX_CHARS];
		if (def_fn != nullptr) {sprintf(fn, "%s", def_fn);}
	}
}

bool load_config_file(const char *fn) {

	string line, config_file;
	ifstream in(fn);
	if (!in.good()) return 0;

	while (std::getline(in, line)) {
		if (line.empty() || line[0] == '#') continue; // comment
		config_file.clear();

		for (auto const &c : line) {
			if (isspace(c) || c == '#') break; // ignore anything after the first string - assume it's a comment
			config_file.push_back(c);
		}
		if (config_file.empty()) continue;
		cout << "Using config file " << config_file << "." << endl;

		if (!load_config(config_file)) {
			cerr << "Error: Failed to open config file " << config_file << " for read" << endl;
			return 0;
		}
	} // end while
	min_eq(NEAR_CLIP, 0.5f*CAMERA_RADIUS); // make sure the player's view can't clip through scene objects
	return 1;
}

int load_top_level_config(const char *def_file) { // defaults.txt

#ifndef _WIN32
	allow_shader_invariants = 0; // Note: I've seen shader invariant errors multiple times on linux but never on Windows, so I guess we disable it by default when not Windows
#endif
	assert(def_file != NULL);
	alloc_if_req(state_file, dstate_file);
	alloc_if_req(mesh_file, dmesh_file);
	alloc_if_req(coll_obj_file, dcoll_obj_file);
	alloc_if_req(ship_def_file, dship_def_file);
	load_config("config_pre.txt"); // defaults
	bool const ret(load_config_file(def_file));
	load_config("config_post.txt"); // overrides (load even if main config file failed)
	apply_grass_scale();
	return ret;
}


void fire_weapon() {

	fire_key = 1;
	if (world_mode == WMODE_INF_TERRAIN) {inf_terrain_fire_weapon();}

	if (world_mode != WMODE_UNIVERSE) {
		assert(sstates != NULL);
		sstates[CAMERA_ID].gamemode_fire_weapon();
	}
}


string get_all_gl_extensions() {

	int n(0);
	glGetIntegerv(GL_NUM_EXTENSIONS, &n);
	string ext;
	for (int i = 0; i < n; i++) {ext += (const char *)glGetStringi(GL_EXTENSIONS, i);}
	return ext;
}

bool has_extension(string const &ext) { // is this always correct?
	return (strstr(get_all_gl_extensions().c_str(), ext.c_str()) != NULL);
}

bool open_file(FILE *&fp, char const *const fn, string const &file_type, char const *const mode) {
	fp = fopen(fn, mode);
	if (fp != nullptr) return 1;
	cerr << "*** Error: Could not open " << file_type << " file '" << fn << "'." << endl;
	return 0;
}

void cfg_err(string const &str, int &error) {
	cerr << "Error reading " << str << " from config file." << endl;
	error = 1;
}


void read_write_lighting_setup(FILE *fp, unsigned ltype, int &error) { // <filename> <write_mode> <light_scale> [<first_ray_weight>]

	assert(ltype < NUM_LIGHTING_TYPES);
	alloc_if_req(lighting_file[ltype], NULL);
	int write_mode(0);
	if (fscanf(fp, "%255s%i%f", lighting_file[ltype], &write_mode, &light_int_scale[ltype]) != 3) {cfg_err("lighting_file command", error);}
	read_float_reset_pos_on_fail(fp, first_ray_weight[ltype]); // ok if fails
	(write_mode ? write_light_files[ltype] : read_light_files[ltype]) = 1;
}


template<typename T> bool kw_to_val_map_t<T>::maybe_set_from_fp(string const &str, FILE *fp) {
	auto it(m.find(str));
	if (it == m.end()) return 0;
	if (!read_type_t(fp, *it->second)) {cfg_err(opt_prefix + str + " keyword", error);}
	return 1;
}
bool kw_to_val_map_float_check_t::map_val_t::check_val() const {
	switch (check_mode) {
	case FP_CHECK_NONE  : return 1;
	case FP_CHECK_POS   : return (*v > 0.0);
	case FP_CHECK_NONNEG: return (*v >= 0.0);
	case FP_CHECK_01    : return (*v >= 0.0 && *v <= 1.0);
	default: assert(0);
	}
	return 1; // never gets here
}
bool kw_to_val_map_float_check_t::maybe_set_from_fp(string const &str, FILE *fp) {
	auto it(m.find(str));
	if (it == m.end()) return 0;
	if (!read_type_t(fp, *it->second.v)) {cfg_err(opt_prefix + str + " keyword", error);}
	if (!it->second.check_val()) {cerr << "Illegal value: " << *it->second.v << "; "; cfg_err(opt_prefix + " " + str + " keyword", error);}
	return 1;
}
template class kw_to_val_map_t<colorRGBA>; // explicit instantiation of this one because it's used for buildings but not here


bool bmp_file_to_binary_array(char const *const fn, unsigned char **&data) {
	if (strlen(fn) > 0) {
		if (!bmp_to_chars(fn, data)) return 0;
		mesh_size_locked = 1;
	}
	return 1;
}


std::string const config_dir("scene_config");

FILE *open_config_file(string const &filename) {

	FILE *fp(fopen(filename.c_str(), "r"));
	if (fp != nullptr) return fp; // found in run dir
	if (open_file(fp, (config_dir + "/" + filename).c_str(), "input configuration file")) return fp; // found in config dir
	return nullptr; // failed
}


int load_config(string const &config_file) {

	FILE *fp(open_config_file(config_file));
	if (fp == nullptr) return 0;
	int error(0);
	char strc[MAX_CHARS] = {0}, md_fname[MAX_CHARS] = {0}, we_fname[MAX_CHARS] = {0}, fw_fname[MAX_CHARS] = {0}, include_fname[MAX_CHARS] = {0};

	// Note: all of these maps bind variable addresses into the config file system by name
	kw_to_val_map_t<bool> kwmb(error);
	kwmb.add("gen_tree_roots", gen_tree_roots);
	kwmb.add("no_smoke_over_mesh", no_smoke_over_mesh);
	kwmb.add("use_waypoints", use_waypoints);
	kwmb.add("use_waypoint_app_spots", use_waypoint_app_spots);
	kwmb.add("group_back_face_cull", group_back_face_cull);
	kwmb.add("inf_terrain_scenery", inf_terrain_scenery);
	kwmb.add("enable_tiled_mesh_ao", enable_tiled_mesh_ao);
	kwmb.add("fast_water_reflect", fast_water_reflect);
	kwmb.add("disable_shader_effects", disable_shader_effects);
	kwmb.add("enable_model3d_tex_comp", enable_model3d_tex_comp);
	kwmb.add("texture_alpha_in_red_comp", texture_alpha_in_red_comp);
	kwmb.add("use_model3d_tex_mipmaps", use_model3d_tex_mipmaps);
	kwmb.add("use_dense_voxels", use_dense_voxels);
	kwmb.add("use_voxel_cobjs", use_voxel_cobjs);
	kwmb.add("mt_cobj_tree_build", mt_cobj_tree_build);
	kwmb.add("global_lighting_update", global_lighting_update);
	kwmb.add("lighting_update_offline", lighting_update_offline);
	kwmb.add("two_sided_lighting", two_sided_lighting);
	kwmb.add("disable_sound", disable_sound);
	kwmb.add("start_maximized", start_maximized);
	kwmb.add("enable_depth_clamp", enable_depth_clamp);
	kwmb.add("detail_normal_map", detail_normal_map);
	kwmb.add("use_core_context", init_core_context);
	kwmb.add("enable_multisample", enable_multisample);
	kwmb.add("dynamic_smap_bias", dynamic_smap_bias);
	kwmb.add("model3d_winding_number_normal", model3d_wn_normal);
	kwmb.add("snow_shadows", snow_shadows);
	kwmb.add("tree_4th_branches", tree_4th_branches);
	kwmb.add("skip_light_vis_test", skip_light_vis_test);
	kwmb.add("model_calc_tan_vect", model_calc_tan_vect);
	kwmb.add("invert_model_nmap_bscale", invert_model_nmap_bscale);
	kwmb.add("enable_dlight_shadows", enable_dlight_shadows);
	kwmb.add("tree_indir_lighting", tree_indir_lighting);
	kwmb.add("only_pine_palm_trees", only_pine_palm_trees);
	kwmb.add("enable_gamma_correction", enable_gamma_correct);
	kwmb.add("use_z_prepass", use_z_prepass);
	kwmb.add("reflect_dodgeballs", reflect_dodgeballs);
	kwmb.add("all_model3d_ref_update", all_model3d_ref_update);
	kwmb.add("store_cobj_accum_lighting_as_blocked", store_cobj_accum_lighting_as_blocked);
	kwmb.add("begin_motion", begin_motion);
	kwmb.add("water_is_lava", water_is_lava);
	kwmb.add("enable_mouse_look", enable_mouse_look);
	kwmb.add("enable_init_shields", enable_init_shields);
	kwmb.add("tt_triplanar_tex", tt_triplanar_tex);
	kwmb.add("enable_model3d_bump_maps", enable_model3d_bump_maps);
	kwmb.add("use_obj_file_bump_grayscale", use_obj_file_bump_grayscale);
	kwmb.add("invert_bump_maps", invert_bump_maps);
	kwmb.add("use_interior_cube_map_refl", use_interior_cube_map_refl);
	kwmb.add("enable_cube_map_bump_maps", enable_cube_map_bump_maps);
	kwmb.add("enable_model3d_custom_mipmaps", enable_model3d_custom_mipmaps);
	kwmb.add("no_store_model_textures_in_memory", no_store_model_textures_in_memory);
	kwmb.add("no_subdiv_model", no_subdiv_model);
	kwmb.add("merge_model_objects", merge_model_objects);
	kwmb.add("use_grass_tess", use_grass_tess);
	kwmb.add("use_instanced_pine_trees", use_instanced_pine_trees);
	kwmb.add("enable_dpart_shadows", enable_dpart_shadows);
	kwmb.add("enable_tt_model_reflect", enable_tt_model_reflect);
	kwmb.add("enable_tt_model_indir", enable_tt_model_indir);
	kwmb.add("auto_calc_tt_model_zvals", auto_calc_tt_model_zvals);
	kwmb.add("disable_tt_water_reflect", disable_tt_water_reflect);
	kwmb.add("use_model_lod_blocks", use_model_lod_blocks);
	kwmb.add("flatten_tt_mesh_under_models", flatten_tt_mesh_under_models);
	kwmb.add("def_texture_compress", def_tex_compress);
	kwmb.add("smileys_chase_player", smileys_chase_player);
	kwmb.add("disable_fire_delay", disable_fire_delay);
	kwmb.add("disable_recoil", disable_recoil);
	kwmb.add("enable_translocator", enable_translocator);
	kwmb.add("enable_grass_fire", enable_grass_fire);
	kwmb.add("tiled_terrain_only", tiled_terrain_only);
	kwmb.add("disable_model_textures", disable_model_textures);
	kwmb.add("start_in_inf_terrain", start_in_inf_terrain);
	kwmb.add("allow_shader_invariants", allow_shader_invariants);
	kwmb.add("unlimited_weapons", config_unlimited_weapons);
	kwmb.add("allow_model3d_quads", allow_model3d_quads);
	kwmb.add("keep_keycards_on_death", keep_keycards_on_death);
	kwmb.add("enable_timing_profiler", enable_timing_profiler);
	kwmb.add("fast_transparent_spheres", fast_transparent_spheres);
	kwmb.add("draw_building_interiors", draw_building_interiors);
	kwmb.add("reverse_3ds_vert_winding_order", reverse_3ds_vert_winding_order);
	kwmb.add("disable_dlights", disable_dlights);
	kwmb.add("enable_hcopter_shadows", enable_hcopter_shadows);
	kwmb.add("pre_load_full_tiled_terrain", pre_load_full_tiled_terrain);
	kwmb.add("disable_blood", disable_blood);
	kwmb.add("enable_model_animations", enable_model_animations);
	kwmb.add("rotate_trees", rotate_trees);
	kwmb.add("invert_model3d_faces", invert_model3d_faces);
	kwmb.add("play_gameplay_alert", play_gameplay_alert);
	kwmb.add("vsync_enabled", vsync_enabled);
	kwmb.add("enable_ground_csm", enable_ground_csm);

	kw_to_val_map_t<int> kwmi(error);
	kwmi.add("verbose", verbose_mode);
	kwmi.add("load_coll_objs", load_coll_objs);
	kwmi.add("glaciate", GLACIATE);
	kwmi.add("dynamic_mesh_scroll", dynamic_mesh_scroll);
	kwmi.add("mesh_seed", mesh_seed);
	kwmi.add("rgen_seed", rgen_seed);
	kwmi.add("universe_only", universe_only);
	kwmi.add("disable_universe", disable_universe);
	kwmi.add("disable_inf_terrain", disable_inf_terrain);
	kwmi.add("left_handed", left_handed);
	kwmi.add("destroy_thresh", destroy_thresh);
	kwmi.add("rand_seed", srand_param);
	kwmi.add("disable_water", INIT_DISABLE_WATER);
	kwmi.add("disable_scenery", DISABLE_SCENERY);
	kwmi.add("read_landscape", read_landscape);
	kwmi.add("read_heightmap", read_heightmap);
	kwmi.add("ground_effects_level", ground_effects_level);
	kwmi.add("tree_coll_level", tree_coll_level);
	kwmi.add("free_for_all", free_for_all);
	kwmi.add("num_dodgeballs", num_dodgeballs);
	kwmi.add("ntrees", num_trees);
	kwmi.add("nsmileys", num_smileys);
	kwmi.add("teams", teams);
	kwmi.add("init_tree_mode", tree_mode);
	kwmi.add("mesh_gen_mode", mesh_gen_mode);
	kwmi.add("mesh_gen_shape", mesh_gen_shape);
	kwmi.add("mesh_freq_filter", mesh_freq_filter);
	kwmi.add("preproc_cube_cobjs", preproc_cube_cobjs);
	kwmi.add("show_waypoints", show_waypoints);
	kwmi.add("init_game_mode", game_mode);
	kwmi.add("init_num_balls", init_num_balls);
	kwmi.add("use_voxel_rocks", use_voxel_rocks); // 0=never, 1=always, 2=only when no vegetation
	kwmi.add("default_anim_id", default_anim_id);

	kw_to_val_map_t<unsigned> kwmu(error);
	kwmu.add("grass_density", grass_density);
	kwmu.add("max_unique_trees", max_unique_trees);
	kwmu.add("shadow_map_sz", shadow_map_sz);
	kwmu.add("max_ray_bounces", MAX_RAY_BOUNCES);
	kwmu.add("num_test_snowflakes", num_snowflakes);
	kwmu.add("hmap_filter_width", hmap_filter_width);
	kwmu.add("erosion_iters", erosion_iters);
	kwmu.add("erosion_iters_tt", erosion_iters_tt);
	kwmu.add("num_dynam_parts", num_dynam_parts);
	kwmu.add("num_birds_per_tile", num_birds_per_tile);
	kwmu.add("num_fish_per_tile", num_fish_per_tile);
	kwmu.add("num_bflies_per_tile", num_bflies_per_tile);
	kwmu.add("max_cube_map_tex_sz", max_cube_map_tex_sz);
	kwmu.add("snow_coverage_resolution", snow_coverage_resolution);
	kwmu.add("dlight_grid_bitshift", DL_GRID_BS);
	kwmu.add("tiled_terrain_gen_heightmap_sz", tiled_terrain_gen_heightmap_sz);
	kwmu.add("game_mode_disable_mask", game_mode_disable_mask);
	kwmu.add("show_map_view_fractal", show_map_view_fractal);

	kw_to_val_map_t<float> kwmf(error);
	kwmf.add("gravity", base_gravity);
	kwmf.add("mesh_height", mesh_height_scale);
	kwmf.add("mesh_scale", mesh_scale);
	kwmf.add("mesh_z_cutoff", mesh_z_cutoff);
	kwmf.add("disabled_mesh_z", disabled_mesh_z);
	kwmf.add("relh_adj_tex", relh_adj_tex);
	kwmf.add("set_czmax", czmax);
	kwmf.add("camera_radius", CAMERA_RADIUS);
	kwmf.add("camera_step_height", C_STEP_HEIGHT);
	kwmf.add("waypoint_sz_thresh", waypoint_sz_thresh);
	kwmf.add("tree_deadness", tree_deadness);
	kwmf.add("tree_dead_prob", tree_dead_prob);
	kwmf.add("sun_rot", sun_rot);
	kwmf.add("moon_rot", moon_rot);
	kwmf.add("sun_theta", sun_theta);
	kwmf.add("moon_theta", moon_theta);
	kwmf.add("cobj_z_bias", cobj_z_bias);
	kwmf.add("indir_vert_offset", indir_vert_offset);
	kwmf.add("self_damage", self_damage);
	kwmf.add("team_damage", team_damage);
	kwmf.add("player_damage", player_damage);
	kwmf.add("smiley_damage", smiley_damage);
	kwmf.add("player_speed", player_speed);
	kwmf.add("smiley_speed", smiley_speed);
	kwmf.add("speed_mult", speed_mult);
	kwmf.add("smiley_accuracy", smiley_acc);
	kwmf.add("crater_size", crater_depth);
	kwmf.add("crater_radius", crater_radius);
	kwmf.add("indir_light_exp", indir_light_exp);
	kwmf.add("snow_random", snow_random);
	kwmf.add("temperature", init_temperature);
	kwmf.add("mesh_start_mag", MESH_START_MAG);
	kwmf.add("mesh_start_freq", MESH_START_FREQ);
	kwmf.add("mesh_mag_mult", MESH_MAG_MULT);
	kwmf.add("mesh_freq_mult", MESH_FREQ_MULT);
	kwmf.add("sm_tree_density", sm_tree_density);
	kwmf.add("tree_density_thresh", tree_density_thresh);
	kwmf.add("tree_slope_thresh", tree_slope_thresh);
	kwmf.add("ocean_wave_height", ocean_wave_height);
	kwmf.add("flower_density", flower_density);
	kwmf.add("model3d_texture_anisotropy", model3d_texture_anisotropy);
	kwmf.add("near_clip_dist", NEAR_CLIP);
	kwmf.add("far_clip_dist", FAR_CLIP);
	kwmf.add("tree_height_scale", tree_height_scale); // applies to trees and small trees
	kwmf.add("sm_tree_scale", sm_tree_scale);
	kwmf.add("model_auto_tc_scale", model_auto_tc_scale);
	kwmf.add("model_triplanar_tc_scale", model_triplanar_tc_scale);
	kwmf.add("smap_thresh_scale", smap_thresh_scale);
	kwmf.add("cloud_height_offset", cloud_height_offset);
	kwmf.add("dodgeball_metalness", dodgeball_metalness);
	kwmf.add("fog_dist_scale", fog_dist_scale);
	kwmf.add("biome_x_offset", biome_x_offset);
	kwmf.add("custom_glaciate_exp", custom_glaciate_exp); // <= 0.0; 0.0 = use default of 3.0
	kwmf.add("tree_type_rand_zone", tree_type_rand_zone); // [0.0, 1.0]
	kwmf.add("universe_ambient_scale", universe_ambient_scale);
	kwmf.add("planet_update_rate", planet_update_rate);
	kwmf.add("jump_height", jump_height);
	kwmf.add("force_czmin", force_czmin);
	kwmf.add("force_czmax", force_czmax);
	kwmf.add("dlight_intensity_scale", dlight_intensity_scale);
	kwmf.add("model_mat_lod_thresh", model_mat_lod_thresh);
	kwmf.add("def_texture_aniso", def_tex_aniso);
	kwmf.add("clouds_per_tile", clouds_per_tile);
	kwmf.add("atmosphere", def_atmosphere);
	kwmf.add("vegetation", def_vegetation);
	kwmf.add("ocean_depth_opacity_mult", ocean_depth_opacity_mult);
	kwmf.add("erode_amount", erode_amount);
	kwmf.add("ambient_scale", ambient_scale);
	kwmf.add("ray_step_size_mult", ray_step_size_mult);
	kwmf.add("system_max_orbit", system_max_orbit);
	kwmf.add("sky_occlude_scale", sky_occlude_scale);
	kwmf.add("mouse_sensitivity", mouse_sensitivity);
	kwmf.add("tt_grass_scale_factor", tt_grass_scale_factor);
	kwmf.add("model_hemi_lighting_scale", model_hemi_lighting_scale);
	kwmf.add("pine_tree_radius_scale", pine_tree_radius_scale);
	kwmf.add("sunlight_brightness", sunlight_brightness);
	kwmf.add("moonlight_brightness", moonlight_brightness);
	kwmf.add("tiled_terrain_fog_density", tt_fog_density); // (0.0, 1.0]
	kwmf.add("mouse_smooth_factor", mouse_smooth_factor); // >= 0.0
	kwmf.add("tree_depth_scale", tree_depth_scale); // >= 0.0

	kwmf.add("hmap_plat_bot",    hmap_params.plat_bot);
	kwmf.add("hmap_plat_height", hmap_params.plat_h);
	kwmf.add("hmap_plat_slope",  hmap_params.plat_s);
	kwmf.add("hmap_plat_max",    hmap_params.plat_max);
	kwmf.add("hmap_crat_height", hmap_params.crat_h);
	kwmf.add("hmap_crat_slope",  hmap_params.crat_s);
	kwmf.add("hmap_crack_lo",    hmap_params.crack_lo);
	kwmf.add("hmap_crack_hi",    hmap_params.crack_hi);
	kwmf.add("hmap_crack_depth", hmap_params.crack_d);
	kwmf.add("hmap_sine_mag",    hmap_params.sine_mag);
	kwmf.add("hmap_sine_freq",   hmap_params.sine_freq);
	kwmf.add("hmap_sine_bias",   hmap_params.sine_bias);
	kwmf.add("hmap_volcano_width",  hmap_params.volcano_width);
	kwmf.add("hmap_volcano_height", hmap_params.volcano_height);

	kw_to_val_map_float_check_t kwmr(error);
	kwmr.add("nleaves_scale",       nleaves_scale,       FP_CHECK_POS);
	kwmr.add("lm_dz_adj",           lm_dz_adj,           FP_CHECK_NONNEG);
	kwmr.add("tree_branch_radius",  branch_radius_scale, FP_CHECK_POS);
	kwmr.add("model3d_alpha_thresh",model3d_alpha_thresh,FP_CHECK_01);
	kwmr.add("snow_depth",          snow_depth,          FP_CHECK_NONNEG);

	kw_to_val_map_t<string> kwms(error);
	kwms.add("cobjs_out_filename", cobjs_out_fn);
	kwms.add("coll_damage_name",   coll_damage_name);
	kwms.add("read_hmap_modmap_filename",  read_hmap_modmap_fn);
	kwms.add("write_hmap_modmap_filename", write_hmap_modmap_fn);
	kwms.add("read_voxel_brush_filename",  read_voxel_brush_fn);
	kwms.add("write_voxel_brush_filename", write_voxel_brush_fn);
	kwms.add("font_texture_atlas_fn", font_texture_atlas_fn);
	kwms.add("sphere_materials_fn", sphere_materials_fn);
	kwms.add("write_heightmap_png", hmap_out_fn);
	kwms.add("skybox_cube_map", skybox_cube_map_name);
	kwms.add("assimp_alpha_exclude_str", assimp_alpha_exclude_str);

	while (read_str(fp, strc)) { // slow but should be OK: these ones require special handling
		string const str(strc);
		if (kwmb.maybe_set_from_fp(str, fp)) continue;
		if (kwmi.maybe_set_from_fp(str, fp)) continue;
		if (kwmu.maybe_set_from_fp(str, fp)) continue;
		if (kwmf.maybe_set_from_fp(str, fp)) continue;
		if (kwmr.maybe_set_from_fp(str, fp)) continue;
		if (kwms.maybe_set_from_fp(str, fp)) continue;

		if (str.size() >= 2 && str[0] == '/' && str[1] == '*') { // start of block comment
			if (!read_block_comment(fp)) {cfg_err("block_comment", error);}
		}
		else if (str[0] == '#') { // comment
			int letter(getc(fp));
			while (letter != '\n' && letter != EOF && letter != 0) letter = getc(fp);
		}
		else if (str == "remap_key") {
			string sfrom, sto;
			if (!read_string(fp, sfrom) || !read_string(fp, sto)) cfg_err("remap_key", error);
			if (!kbd_remap.parse_remap_command(sfrom, sto)) cfg_err("remap_key", error);
		}
		else if (str == "voxel") { // voxel option
			if (!parse_voxel_option(fp)) cfg_err("voxel option", error);
		}
		else if (str == "buildings") { // buildings option
			if (!parse_buildings_option(fp)) cfg_err("buildings option", error);
		}
		else if (str == "city") { // city options
			if (!parse_city_option(fp)) cfg_err("city option", error);
		}
		else if (str == "sphere_gen") { // sphere_gen options
			if (!parse_sphere_gen_option(fp)) cfg_err("sphere_gen option", error);
		}
		else if (str == "include") {
			if (!read_str(fp, include_fname)) cfg_err("include", error);
			if (!load_config(include_fname )) cfg_err("nested include file", error);
		}
		else if (str == "grass_size") {
			if (!read_pos_float(fp, grass_length) || !read_pos_float(fp, grass_width)) {cfg_err("grass size", error);}
		}
		else if (str == "force_tree_class") {
			if (!read_int(fp, force_tree_class) || force_tree_class >= NUM_TREE_CLASSES) cfg_err("force_tree_class", error);
		}
		else if (str == "tree_lod_scale") {
			for (unsigned i = 0; i < 4; ++i) {
				if (!read_non_neg_float(fp, tree_lod_scales[i])) cfg_err("tree_lod_scale", error);
			}
			if (tree_lod_scales[0] < tree_lod_scales[1] || tree_lod_scales[2] < tree_lod_scales[3]) {cfg_err("tree_lod_scale values", error);}
		}
		else if (str == "num_items") { // HEALTH, SHIELD, POWERUP, WEAPON, AMMO
			for (unsigned n = 0; n < sizeof(init_item_counts)/sizeof(unsigned); ++n) {
				if (!read_uint(fp, init_item_counts[n])) {cfg_err("number of items", error); break;}
			}
		}
		else if (str == "window_width") {
			if (!read_int(fp, window_width ) || window_width  < 1) cfg_err("window_width command", error);
		}
		else if (str == "window_height") {
			if (!read_int(fp, window_height) || window_height < 1) cfg_err("window_height command", error);
		}
		else if (str == "mesh_size") {
			if (fscanf(fp, "%i%i%i", &MESH_X_SIZE, &MESH_Y_SIZE, &MESH_Z_SIZE) != 3) cfg_err("mesh size command", error);

			if (mesh_size_locked) {
				cerr << "Error: mesh_size command cannot be called after loading a file that depends on the size" << endl;
				error = 1;
			}
		}
		else if (str == "scene_size") {
			if (fscanf(fp, "%f%f%f", &X_SCENE_SIZE, &Y_SCENE_SIZE, &Z_SCENE_SIZE) != 3) cfg_err("scene size command", error);
		}
		else if (str == "load_hmv") {
			if (!read_int(fp, load_hmv)) cfg_err("load_hmv command", error);
			if (fscanf(fp, "%f%f%f%f", &hmv_pos.x, &hmv_pos.y, &hmv_pos.z, &hmv_scale) != 4) cfg_err("load_hmv command", error);
		}
		else if (str == "water_h_off") { // abs [rel]
			if (!read_float(fp, water_h_off)) cfg_err("water_h_off command", error);
			water_h_off_rel = 0.0;
			read_float(fp, water_h_off_rel); // optional
		}
		else if (str == "wind_velocity") {
			if (!read_vector(fp, wind)) cfg_err("wind_velocity command", error);
		}
		else if (str == "camera_height") {
			if (!read_double(fp, camera_zh)) cfg_err("camera_height command", error);
		}
		else if (str == "player_start") {
			player_custom_start_pos = 1;
			if (!read_vector(fp, surface_pos)) cfg_err("player_start command", error);
		}
		else if (str == "cube_map_center") {
			if (!read_vector(fp, cube_map_center)) cfg_err("cube_map_center command", error);
			if (cube_map_center == all_zeros) {cube_map_center.x += TOLERANCE;} // since all_zeros is a special flag, make sure the user-specified value is different
		}
		else if (str == "tree_size") {
			float tree_size(1.0);
			if (!read_pos_float(fp, tree_size)) cfg_err("tree size command", error);
			tree_scale = 1.0/tree_size;
		}
		else if (str == "bush_probability") {
			for (unsigned i = 0; i < NUM_TREE_TYPES; ++i) { // read a list of floating-point numbers (could allow a partial set to be read)
				if (!read_zero_one_float(fp, tree_types[i].bush_prob)) {cfg_err("bush_probability command", error); break;}
			}
		}
		else if (str == "leaf_color") {
			if (fscanf(fp, "%f%f%f%f%f", &leaf_base_color.R, &leaf_base_color.G, &leaf_base_color.B, &leaf_color_coherence, &tree_color_coherence) != 5) {
				cfg_err("leaf_color command", error);
			}
		}
		else if (str == "flower_color") {
			if (fscanf(fp, "%f%f%f", &flower_color.R, &flower_color.G, &flower_color.B) != 3) {cfg_err("flower_color command", error);}
			flower_color.A = 1.0;
		}
		else if (str == "sunlight_color") {
			if (fscanf(fp, "%f%f%f", &sunlight_color.R, &sunlight_color.G, &sunlight_color.B) != 3) {cfg_err("sunlight_color command", error);}
		}
		else if (str == "sunlight_intensity") {
			float intensity(1.0);
			if (!read_float(fp, intensity)) {cfg_err("sunlight_intensity command", error);}
			sunlight_color *= intensity;
		}
		else if (str == "floating_light_params") { // rmin, rmax, vmin, vmax, imin, imax
			if (fscanf(fp, "%f%f%f%f%f%f", &dp_params.rmin, &dp_params.rmax, &dp_params.vmin, &dp_params.vmax, &dp_params.imin, &dp_params.imax) != 6) {
				cfg_err("floating_light_params command", error);
			}
		}
		else if (str == "floating_light_range") { // x1 x2 y1 y2 z1 z2
			if (fscanf(fp, "%f%f%f%f%f%f", &dp_params.sdist[0].x, &dp_params.sdist[1].x, &dp_params.sdist[0].y, &dp_params.sdist[1].y, &dp_params.sdist[0].z, &dp_params.sdist[1].z) != 6) {
				cfg_err("floating_light_range command", error);
			}
		}
		else if (str == "toggle_mesh_enabled") {display_mode ^= 0x01;}
		else if (str == "toggle_reflections" ) {display_mode ^= 0x10;}
		else if (str == "player_name") {
			if (!read_str(fp, player_name)) cfg_err("player name", error);
		}
		else if (str == "create_voxel_landscape") {
			if (!read_uint(fp, create_voxel_landscape) || create_voxel_landscape > 2) cfg_err("create_voxel_landscape command", error);
		}
		else if (str == "team_start") {
			bbox bb;
			if (fscanf(fp, "%i%f%f%f%f", &bb.index, &bb.x1, &bb.y1, &bb.x2, &bb.y2) != 5) cfg_err("team start command", error);
			if (bb.index < 0 || bb.index >= teams) {cout << "Error: Illegal team specified in team_start command: " << bb.index << " (" << teams << " teams)." << endl; error = 1;}
			else team_starts.push_back(bb);
		}
		else if (str == "vertex_optimize_flags") {
			for (unsigned i = 0; i < 3; ++i) { // {enable, full_opt, verbose}
				if (!read_bool(fp, vert_opt_flags[i])) cfg_err("vertex_optimize_flags command", error);
			}
		}
		else if (str == "coll_obj_file") {
			if (!read_str(fp, coll_obj_file)) cfg_err("coll_obj_file command", error);
		}
		else if (str == "state_file") {
			if (!read_str(fp, state_file)) cfg_err("state_file command", error);
		}
		else if (str == "mesh_file") { // only the first parameter is required
			float rmz(0.0);
			if (fscanf(fp, "%255s%f%f%i", mesh_file, &mesh_file_scale, &mesh_file_tz, &do_read_mesh) < 1) cfg_err("mesh_file command", error);
			if (read_float(fp, rmz)) {read_mesh_zmm = rmz;}
		}
		else if (str == "mh_filename") { // only the first parameter is required
			alloc_if_req(mh_filename, NULL);
			if (fscanf(fp, "%255s%f%f%i", mh_filename, &mesh_file_scale, &mesh_file_tz, &invert_mh_image) < 1) cfg_err("mh_filename command", error);
		}
		else if (str == "mh_filename_tiled_terrain") {
			alloc_if_req(mh_filename_tt, NULL);
			if (fscanf(fp, "%255s", mh_filename_tt) != 1) cfg_err("mh_filename_tiled_terrain command", error);
		}
		else if (str == "mesh_diffuse_tex_fn") {
			alloc_if_req(mesh_diffuse_tex_fn, NULL);
			if (fscanf(fp, "%255s", mesh_diffuse_tex_fn) != 1) cfg_err("mesh_diffuse_tex_fn command", error);
			read_bool(fp, mesh_difuse_tex_comp); // okay if fails
		}
		else if (str == "default_ground_tex") {
			if (!read_str(fp, strc)) cfg_err("default_ground_tex", error);
			default_ground_tex = get_texture_by_name(string(strc));
		}
		else if (str == "mesh_detail_tex") {
			if (!read_str(fp, strc)) cfg_err("mesh_detail_tex", error);
			mesh_detail_tex = get_texture_by_name(string(strc));
		}
		else if (str == "skybox_tex") {
			if (!read_str(fp, strc)) cfg_err("skybox_tex", error);
			string const texture_name = string(strc);
			if (!check_texture_file_exists(texture_name)) {std::cerr << "Error: Failed to load skybox texture '" << texture_name << "'; Disabling skybox" << endl;}
			else {skybox_tid = get_texture_by_name(texture_name, 0, 0, 0);} // clamp
		}
		else if (str == "ship_def_file") {
			if (!read_str(fp, ship_def_file)) cfg_err("ship_def_file command", error);
		}
		else if (str == "smoke_bounds") {
			cube_t sb;
			for (unsigned d = 0; d < 6; ++d) { // x1 x2 y1 y2 z1 z2
				if (!read_float(fp, sb.d[d>>1][d&1])) cfg_err("smoke_bounds command", error);
			}
			smoke_bounds.push_back(sb);
		}
		else if (str == "reflect_plane_z") {
			cube_t cube;
			if (read_cube(fp, geom_xform_t(), cube) != 6) cfg_err("reflect_plane_z command", error);
			reflect_planes.add(cube);
		}
		// lighting
		else if (str == "lighting_file_sky") {
			read_write_lighting_setup(fp, LIGHTING_SKY, error); // <filename> <write_mode> <scale>
		}
		else if (str == "lighting_file_global") {
			read_write_lighting_setup(fp, LIGHTING_GLOBAL, error); // <filename> <write_mode> <scale> [<first ray weight>]
		}
		else if (str == "lighting_file_local") {
			read_write_lighting_setup(fp, LIGHTING_LOCAL, error); // <filename> <write_mode> <scale>
		}
		else if (str == "lighting_file_cobj") {
			read_write_lighting_setup(fp, LIGHTING_COBJ_ACCUM, error); // <filename> <write_mode> <scale>
		}
		else if (str == "num_light_rays") { // GLOBAL_RAYS and DYNAMIC_RAYS are optional
			if (fscanf(fp, "%u%u%u%u%u", &NPTS, &NRAYS, &LOCAL_RAYS, &GLOBAL_RAYS, &DYNAMIC_RAYS) < 3) cfg_err("num_light_rays command", error);
		}
		else if (str == "num_threads") {
			if (!read_uint(fp, NUM_THREADS) || NUM_THREADS > 100) cfg_err("num_threads", error);
			
			if (NUM_THREADS == 0) { // auto special case
				unsigned const num_hw_threads(std::thread::hardware_concurrency()); // includes hyperthreading
				NUM_THREADS = ((num_hw_threads > 0) ? num_hw_threads : 4); // default to 4 if the call returns 0
			}
		}
		else if (str == "ambient_lighting_scale") {
			if (fscanf(fp, "%f%f%f", &ambient_lighting_scale.R, &ambient_lighting_scale.G, &ambient_lighting_scale.B) != 3) cfg_err("ambient_lighting_scale command", error);
		}
		else if (str == "mesh_color_scale") {
			if (fscanf(fp, "%f%f%f", &mesh_color_scale.R, &mesh_color_scale.G, &mesh_color_scale.B) != 3) cfg_err("mesh_color_scale command", error);
		}
		// snow
		else if (str == "snow_file") {
			alloc_if_req(snow_file, NULL);
			int write_mode(0);
			if (fscanf(fp, "%255s%i", snow_file, &write_mode) != 2) cfg_err("snow_file command", error);
			(write_mode ? write_snow_file : read_snow_file) = 1;
		}
		// image files
		else if (str == "mesh_draw_bmp") {
			if (!read_str(fp, md_fname)) cfg_err("mesh_draw_bmp command", error);
		}
		else if (str == "water_enabled_bmp") {
			if (!read_str(fp, we_fname)) cfg_err("water_enabled_bmp command", error);
		}
		else if (str == "flower_weight_bmp") {
			if (!read_str(fp, fw_fname)) cfg_err("flower_weight_bmp command", error);
		}
		else if (str == "end") {
			break;
		}
		else {
			cout << "Unrecognized keyword in input file: " << str << endl;
			error = 1;
		}
		if (error) {cout << "Parse error in config file." << endl; break;}
	} // while read
	if (universe_only && disable_universe) {cout << "Error: universe_only and disable_universe are mutually exclusive" << endl; error = 1;}
	if (mh_filename_tt != nullptr && tiled_terrain_gen_heightmap_sz > 0) {cout << "Error: can't specify both mh_filename_tiled_terrain and tiled_terrain_gen_heightmap_sz" << endl; error = 1;}
	checked_fclose(fp);
	temperature    = init_temperature;
	num_dodgeballs = max(num_dodgeballs, 1); // have to have at least 1
	num_trees      = max(num_trees,      0);
	num_smileys    = max(num_smileys,    0);
	teams          = max(teams,          1);
	tree_mode      = tree_mode % 4;
	use_core_context = init_core_context;
	if (universe_only) {world_mode = WMODE_UNIVERSE;}
	if (tiled_terrain_only || start_in_inf_terrain) {world_mode = WMODE_INF_TERRAIN;}
	//if (read_heightmap && dynamic_mesh_scroll) cout << "Warning: read_heightmap and dynamic_mesh_scroll are currently incompatible options as the heightmap does not scroll." << endl;
	DISABLE_WATER = INIT_DISABLE_WATER;
	XY_MULT_SIZE  = MESH_X_SIZE*MESH_Y_SIZE; // for bmp_to_chars() allocation
	precip_dist_scale = CAMERA_RADIUS/DEF_CAMERA_RADIUS;
	if (!bmp_file_to_binary_array(md_fname, mesh_draw    )) {error = 1;}
	if (!bmp_file_to_binary_array(we_fname, water_enabled)) {error = 1;}
	if (!bmp_file_to_binary_array(fw_fname, flower_weight)) {error = 1;}
	if (!error && !sphere_materials_fn.empty()) {error = !read_sphere_materials_file(sphere_materials_fn);}
	if (error) exit(1);
	return 1;
}


void progress() {cout << "."; cout.flush();}


#ifdef _WIN32
void APIENTRY
#else
void
#endif
openglCallbackFunction(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {

	if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return; // don't spam stdout with notifications
	if (type == GL_DEBUG_TYPE_PERFORMANCE && severity == GL_DEBUG_SEVERITY_MEDIUM && id == 131218) return; // perf warning: shader recompile based on state change, skip
	cout << "### OpenGL message: " << message << endl << "type: ";
	switch (type) {
	case GL_DEBUG_TYPE_ERROR: cout << "ERROR"; break;
	case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: cout << "DEPRECATED_BEHAVIOR"; break;
	case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: cout << "UNDEFINED_BEHAVIOR"; break;
	case GL_DEBUG_TYPE_PORTABILITY: cout << "PORTABILITY"; break;
	case GL_DEBUG_TYPE_PERFORMANCE: cout << "PERFORMANCE"; break;
	case GL_DEBUG_TYPE_OTHER: cout << "OTHER"; break;
	}
	cout << " id: " << id << " severity: ";
	switch (severity) {
	case GL_DEBUG_SEVERITY_LOW: cout << "LOW"; break;
	case GL_DEBUG_SEVERITY_MEDIUM: cout << "MEDIUM"; break;
	case GL_DEBUG_SEVERITY_HIGH: cout << "HIGH"; break;
	case GL_DEBUG_SEVERITY_NOTIFICATION: cout << "NOTIFICATION"; break;
	default: cout << hex << severity << dec << " ";
	}
	cout << endl;
	//assert(type != GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR); // uncomment to break/exit on this type of error
	//assert(severity != GL_DEBUG_SEVERITY_HIGH); // uncomment to break/exit on this type of error
}

void init_debug_callback() {

#if _DEBUG
	if (glDebugMessageCallback) {
		cout << "Register OpenGL debug callback " << endl;
		glEnable(GL_DEBUG_OUTPUT); // is this needed?
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		glDebugMessageCallback(openglCallbackFunction, nullptr);
		GLuint unusedIds = 0;
		glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, &unusedIds, true);
	}
	else {cout << "glDebugMessageCallback not available" << endl;}
#endif
}


int main(int argc, char** argv) {

	cout << "Starting 3DWorld" << endl;
	//HeapSetInformation(NULL, HeapEnableTerminationOnCorruption, NULL, 0);
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	if (argc == 2) {read_ueventlist(argv[1]);}
	int rs(1);
	if      (srand_param == 1) {rs = GET_TIME_MS();}
	else if (srand_param != 0) {rs = srand_param;}
	add_uevent_srand(rs);
	create_sin_table();
	set_scene_constants();
	load_texture_names(); // needs to be before config file load
	load_top_level_config(defaults_file);
	gen_gauss_rand_arr(); // after reading seed from config file
	cout << "Loading."; cout.flush();
	
 	// Initialize GLUT
	progress();
 	glutInit(&argc, argv);
	progress();
 	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
	glutInitWindowSize(window_width, window_height);

	if (init_core_context) {
		glutInitContextVersion(4, 5);
		glutInitContextFlags(GLUT_CORE_PROFILE
#if _DEBUG
			| GLUT_DEBUG
#endif
		);
	}
	if (enable_timing_profiler) {toggle_timing_profiler();} // enable profiler logging on init without using the 'u' key
	progress();
	glutCreateWindow("3D World");
	progress();
	init_openal(argc, argv);
	progress();
	init_glew();
	progress();
	init_window();
	check_gl_error(7770);
	if (init_core_context) {init_debug_callback();}
	//glEnable(GL_FRAMEBUFFER_SRGB);
	cout << ".GL Initialized." << endl;
	uevent_advance_frame();
	--frame_counter;

	if (0) { // reversed Z-buffer; OpenGL 4.5 only; appears to break shadows
		// see https://www.wedesoft.de/software/2021/09/20/reversed-z-rendering/
		glClipControl(GL_LOWER_LEFT, GL_ZERO_TO_ONE);
	}
	check_gl_error(7771);
	load_textures();
	load_flare_textures(); // Sun Flare
	check_gl_error(7772);
	setup_shaders();
	check_gl_error(7773);
	//cout << "Extensions: " << get_all_gl_extensions() << endl;

	if (!universe_only) { // universe mode should be able to do without these initializations
		reset_planet_defaults(); // set atmosphere and vegetation
		init_objects();
		alloc_matrices();
		t_trees.resize(num_trees);
		init_models();
		init_terrain_mesh();
		init_lights();
		check_gl_error(7774);
		gen_scene(1, (world_mode == WMODE_GROUND), 0, 0, 0);
		check_gl_error(7775);
		gen_snow_coverage();
		if (enable_grass_fire) {init_ground_fire();}
		create_object_groups();
		init_game_state();
		check_gl_error(7776);

		if (game_mode && world_mode == WMODE_GROUND) {
			gamemode_rand_appear();
			camera_mode = 1; // on the ground
		}
		get_landscape_texture_color(0, 0); // hack to force creation of the cached_ls_colors vector in the master thread (before build_lightmap())
		build_lightmap(1);
	}
	check_gl_error(7777);
	glutMainLoop(); // Switch to main loop
	quit_3dworld(); // never actually gets here
    return 0;
}

