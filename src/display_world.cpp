// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "mesh.h"
#include "main.h"
#include "timetest.h"
#include "physics_objects.h"
#include "profiler.h"
#include "model3d.h"
#include <fstream>


/* lights used:
0: Sun, Universe Star
1: Moon, Universe Ambient (universe mode)
2: Lightning, Uobjs
3: Universe Ambient (combined_gu), Uobjs
4: Moon and Sky (local), Uobjs
5: Uobjs
6: Ship Engines, Uobjs
7: Ship Engines */


//#define USE_GPU_TIMER
float const DSCALE             = 0.45;
float const ASCALE             = 0.35;
bool  const TIMETEST           = (GLOBAL_TIMETEST || 0);
bool  const DETERMINISTIC_TIME = 0;
int   const KEEP_MESH          = 1;
float const REL_SCROLL_DIST    = 0.4;
float const LIGHT_W_VAL        = 0.0; // 0 = directional/light at infinity, 1 = point source
float const C_RADIUS0          = 0.01;
float const CR_SCALE           = 0.1;
float const FOG_COLOR_ATTEN    = 0.75;


bool mesh_invalidated(1), fog_enabled(0), tt_fire_button_down(0), in_loading_screen(0), no_tt_footsteps(0);
int iticks(0), time0(0), scrolling(0), dx_scroll(0), dy_scroll(0);
unsigned enabled_lights(0), cur_display_iter(0); // 8 bit flags for enabled_lights
float fticks(0.0), tstep(0.0), camera_shake(0.0), cur_fog_end(1.0), far_clip_ratio(1.0);
double tfticks(0.0), sim_ticks(0.0);
upos_point_type cur_origin(all_zeros);
colorRGBA bkg_color(LT_BLUE), cur_fog_color(GRAY), base_cloud_color(WHITE), base_sky_color(BACKGROUND_DAY), sunlight_color(SUN_LT_C);
string lighting_update_text;


extern bool combined_gu, have_sun, clear_landscape_vbo, show_lightning, spraypaint_mode, enable_depth_clamp, enable_multisample, water_is_lava;
extern bool user_action_key, flashlight_on, enable_clip_plane_z, begin_motion, config_unlimited_weapons, start_maximized, show_bldg_pickup_crosshair;
extern bool can_do_building_action, enable_tt_model_indir, pre_load_full_tiled_terrain, player_custom_start_pos, player_in_tunnel;
extern unsigned inf_terrain_fire_mode, reflection_tid;
extern int auto_time_adv, camera_flight, reset_timing, run_forward, window_width, window_height, voxel_editing, UNLIMITED_WEAPONS, player_in_basement;
extern int advanced, b2down, dynamic_mesh_scroll, spectate, animate2, used_objs, disable_inf_terrain, DISABLE_WATER, can_pickup_bldg_obj, player_in_water;
extern float TIMESTEP, NEAR_CLIP, FAR_CLIP, cloud_cover, univ_sun_rad, atmosphere, vegetation, zmin, zbottom, ztop, ocean_wave_height, brightness;
extern float def_atmosphere, def_vegetation, clip_plane_z, ambient_scale, sunlight_brightness, moonlight_brightness;
extern point mesh_origin, surface_pos, univ_sun_pos, orig_cdir, sun_pos, moon_pos, debug_event_pos;
extern vector3d total_wind;
extern colorRGBA sun_color;
extern water_params_t water_params;
extern lightning_t l_strike;
extern vector<camera_filter> cfilters;
extern reflective_cobjs_t reflective_cobjs;

void check_xy_offsets();
void post_window_redisplay();
void display_universe();
void display_inf_terrain();
void update_temperature(bool verbose);
void update_sound_loops();
bool indir_lighting_updated();
point get_universe_display_camera_pos();
colorRGBA get_inf_terrain_mod_color();
void run_postproc_effects();
void play_camera_footstep_sound();
void draw_voxel_edit_volume();
void play_switch_weapon_sound();
void toggle_fullscreen();
void calc_cur_ambient_diffuse();
void print_texture_memory_usage();
void show_debug_event_pos(point const &pos);

vector3d calc_camera_direction();
void draw_player_model(point const &pos, vector3d const &dir, int time);
void building_gameplay_next_frame(); // from building_interact.cc
void create_building_reflections();
void follow_city_actor();


void glClearColor_rgba(const colorRGBA &color) {
	glClearColor(color.R, color.G, color.B, color.A);
}


void calc_moon_atten(colorRGBA &ambient, colorRGBA &diffuse, float mlf) {

	if (mlf < 0.6) {diffuse *= ((mlf < 0.5) ? 0.0 : 10.0*(mlf - 0.5));}
	ambient *= 0.5;
	diffuse *= 0.5;
}


void set_standard_viewport() {
	glViewport(0, 0, window_width, window_height);
}

void set_player_pdu(vector3d const &rv1, vector3d const &rv2) {

	vector3d cview_dir_n(cview_dir), upv(up_vector);
	rotate_vector3d_by_vr(rv1, rv2, cview_dir_n); //if (rv1 != rv2)?
	rotate_vector3d_by_vr(rv1, rv2, upv);
	set_player_up(upv);
	set_player_dir(cview_dir_n.get_norm());
	set_univ_pdu();
}

void do_look_at() {

	point eye, center;

	if (world_mode == WMODE_UNIVERSE) { // special universe pre-translation to reduce floating-point magnitudes and errors
		up_vector     = get_player_up();
		cur_origin    = get_player_pos();
		camera_origin = all_zeros;
		//camera_origin = get_player_pos(); // combined_gu?
		eye    = camera_origin;
		center = eye + get_player_dir();
		set_univ_pdu();
		up_vector.x += 0.05*camera_shake;
		up_vector.normalize();
	}
	else {
		set_camera_pdu();
		assert(cview_radius > TOLERANCE);
		assert(cview_dir != zero_vector);
		cur_origin = all_zeros;
		center     = camera_origin;
		center.x  += 0.05*CAMERA_RADIUS*camera_shake;
		eye        = center - cview_dir*cview_radius;
	}
	camera_shake = -pow(0.95f, fticks)*camera_shake;
	if (fabs(camera_shake) < 0.1) {camera_shake = 0.0;}
	assert(!dist_less_than(eye, center, TOLERANCE));
	assert(up_vector != zero_vector);
	fgLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up_vector.x, up_vector.y, up_vector.z);
}


void apply_camera_offsets(point const &camera) { // for TT mode

	int const dx(int(camera.x*DX_VAL_INV) - xoff);
	int const dy(int(camera.y*DY_VAL_INV) - yoff);

	if (dx != 0 || dy != 0) {
		xoff  += dx;
		yoff  += dy;
		xoff2 += dx;
		yoff2 += dy;
	}
}


point get_sun_pos() {

	point pos(sun_pos);
	// Note: moving the sun and moon along with the camera makes the sun seem to be infintely far away and looks more real;
	// however, it also makes the shadows misalign slightly, especially for large scenes;
	// we can't have the shadows move with the camera either because then we would have to recompute them every time the player moves
	// Note: if this line is commented out, we must also change inf terrain mode sun drawing
	if (camera_mode == 1) {pos += surface_pos;}
	if (camera_mode != 1 && combined_gu) {pos += get_camera_pos();} // universe is always centered around the camera
	return pos;
}


point get_moon_pos() {

	point pos(moon_pos);
	if (camera_mode == 1) {pos += surface_pos;}
	return pos;
}


colorRGBA get_bkg_color(point const &p1, vector3d const &v12) { // optimize?

	colorRGBA color(bkg_color);
	// use textures[GRADIENT_TEX]?
	
	if (light_factor > 0.4) {
		point const spos(get_sun_pos());
		float const outer_radius(4.0*sun_radius);
		float rad, t, dist;

		if (line_intersect_sphere(p1, v12, spos, outer_radius, rad, dist, t)) {
			if (rad < sun_radius) {color = sun_color;}
			else {
				float blend_const(1.0 - (rad - sun_radius)/(outer_radius - sun_radius)); // 1.0 = full sun, 0.0 = no sun
				blend_color(color, sun_color, color, blend_const*blend_const, 0);
			}
		}
	}
	if (combined_gu) return color;
	if (light_factor < 0.6 && line_intersect_sphere(p1, v12, get_moon_pos(), moon_radius)) {color = texture_color(MOON_TEX);}
	return color;
}


void draw_stuff(int draw_uw, int timer1, int reflection_pass=0) {

	if (draw_uw) {
		draw_bubbles();
		if (underwater) {draw_underwater_particles(min(zmin, zbottom));}
	}
	else { // camera above water
		draw_splashes();
		draw_snow();
		draw_trees(0, (reflection_pass != 0));
		render_voxel_data(0);
		check_gl_error(20);
		if (TIMETEST) PRINT_TIME("O");
		render_models(0, reflection_pass);
		check_gl_error(21);
		if (TIMETEST) PRINT_TIME("P");
		draw_jump_pads();
		draw_teleporters();
		if (show_lightning) {l_strike.gen();} // after the water?

		if (!underwater) {
			maybe_draw_rainbow();

			if (reflection_pass != 1) { // don't draw precip in planar reflections
				draw_local_precipitation(reflection_pass != 0);
				check_gl_error(23);
				if (TIMETEST) PRINT_TIME("R");
			}
		}
		//draw_spotlight_cones();
		draw_coll_surfaces(1, reflection_pass); // transparent
		check_gl_error(24);
		if (TIMETEST) PRINT_TIME("S");
		draw_cracks_and_decals();
		if (TIMETEST) PRINT_TIME("S2");
		draw_transparent_object_groups(reflection_pass);
		draw_voxel_edit_volume();
		check_gl_error(25);
	}
}


void log_location(point const &pos) {

	static std::ofstream out;
	static bool inited(0);

	if (!inited) {
		out.open("positions.log.txt");
		inited = 1;
	}
	assert(out.good());
	out << pos.x << " " << pos.y << " " << pos.z << endl;
}


void draw_frame_rate(float framerate) {

	if (show_framerate) {
		point const camera((world_mode == WMODE_UNIVERSE) ? get_universe_display_camera_pos() : get_camera_pos());
		cout << "FPS: " << framerate << "  loc: (" << camera.str() << ") @ frame " << frame_counter << endl;
		log_location(camera);
		show_framerate = 0;
	}
	if (display_framerate) {
		static int fr_counter(0);
		static float fr2(0.0);
		static double last_tfticks(0.0);
		if (tfticks > last_tfticks + 0.1*TICKS_PER_SECOND) {fr2 = framerate; last_tfticks = tfticks;} // update every 100ms
		draw_framerate(fr2);
		++fr_counter;
	}
}

class framerate_tracker_t {
	int last_report_frame=0;
	float fr_average=0.0;
	high_resolution_clock::time_point timer_a, global_time, last_report_time, timer_b;
	high_resolution_clock clock;
public:
	float get_framerate() {
		timer_b = clock.now();

		if (timer_b > timer_a) { // skip zero time frames
			float const framerate(1.0f/get_delta_secs(timer_b, timer_a));
			timer_a = timer_b;
			//return framerate;
			float const NUM_AVG = 5; // average over several frames
			fr_average = ((fr_average == 0.0) ? framerate : ((NUM_AVG - 1)*fr_average + framerate)/NUM_AVG);
		}
		return fr_average;
	}
	float show_cur_framerate() {
		float const framerate(get_framerate());

		if (show_framerate == 2) {
			cout << used_objs << " objects, time = " << 1000.0*get_delta_secs(timer_b, global_time) << endl;
			cout << "Elapsed frames = " << (frame_counter - last_report_frame) << ", elapsed time = " << 1000.0*get_delta_secs(timer_b, last_report_time)
				<< ", avg framerate = " << ((float)frame_counter - (float)last_report_frame)/get_delta_secs(timer_b, last_report_time) << endl;
			last_report_frame = frame_counter;
			last_report_time  = timer_b;
			show_framerate    = 0;
		}
		if (show_framerate == 1) {
			timer_a = clock.now();
			show_framerate = 2;
		}
		return framerate;
	}
};
framerate_tracker_t framerate_tracker;


// Nvidia specific defines
#define GPU_MEMORY_INFO_DEDICATED_VIDMEM_NVX          0x9047
#define GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX    0x9048
#define GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX  0x9049
#define GPU_MEMORY_INFO_EVICTION_COUNT_NVX            0x904A
#define GPU_MEMORY_INFO_EVICTED_MEMORY_NVX            0x904B
// AMD/ATI specific defines
#define VBO_FREE_MEMORY_ATI                           0x87FB
#define TEXTURE_FREE_MEMORY_ATI                       0x87FC
#define RENDERBUFFER_FREE_MEMORY_ATI                  0x87FD

void show_gpu_mem_info() {
	check_gl_error(10111); // in case there was an incoming error
	int ded_vmem(0), tot_vmem(0), avail_vmem(0), vbo_free_mem(0), texture_free_mem(0), rbuf_free_mem(0);
	glGetIntegerv(GPU_MEMORY_INFO_DEDICATED_VIDMEM_NVX,         &ded_vmem);
	glGetIntegerv(GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX,   &tot_vmem);
	glGetIntegerv(GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX, &avail_vmem);
	glGetIntegerv(VBO_FREE_MEMORY_ATI,                          &vbo_free_mem);
	glGetIntegerv(TEXTURE_FREE_MEMORY_ATI,                      &texture_free_mem);
	glGetIntegerv(RENDERBUFFER_FREE_MEMORY_ATI,                 &rbuf_free_mem);
	glGetError(); // clear the error state
	if (ded_vmem    ) {cout << TXT(ded_vmem) << TXT(tot_vmem) << TXT(avail_vmem) << endl;} // Nvidia
	if (vbo_free_mem) {cout << TXT(vbo_free_mem) << TXT(texture_free_mem) << TXT(rbuf_free_mem) << endl;} // ATI
}


void final_draw(float framerate) {

	fog_enabled = 0;
	run_postproc_effects();
	check_zoom(); // also resets MVM to identity
	draw_inventory(); // drawn last, on top of everything else
	draw_camera_filters(cfilters);
	draw_frame_rate(framerate);
	show_other_messages();
	user_action_key = 0;
	//show_gpu_mem_info(); // TESTING
}


void swap_buffers_and_redraw() {
	glutSwapBuffers();
	if (animate) {post_window_redisplay();} // before glutSwapBuffers()?
}


float get_lf_scale(float lf)   {return CLIP_TO_01(5.0f*(lf - 0.4f));}
float get_light_factor_scale() {return get_lf_scale(light_factor);}
float get_moon_light_factor()  {return fabs(moon_rot/PI - 1.0);}


float get_star_alpha(bool obscured_by_clouds) {

	float star_alpha(obscured_by_clouds ? 0.0f : (1.0f - get_light_factor_scale()));
	if (star_alpha >= 1.0) {return 1.0;}
	//if (world_mode != WMODE_INF_TERRAIN) {return star_alpha;}
	float const dist_above_clouds(get_camera_pos().z - get_tt_cloud_level());
	if (dist_above_clouds > 0.0) {star_alpha = min(1.0, (star_alpha + 0.05*dist_above_clouds));}
	star_alpha = atmosphere*star_alpha + (1.0 - atmosphere); // star alpha increases as atmosphere decreases
	return star_alpha;
}

colorRGBA attenuate_sun_color(colorRGBA const &c) {

	colorRGBA sc;
	sc.A = 1.0;
	UNROLL_3X(sc[i_] = c[i_]*sunlight_color[i_]/SUN_LT_C[i_];) // scale based on ratio of sunlight color to default value
	sc.normalize_to_max_comp();
	return sc;
}

void calc_bkg_color() {

	float const star_alpha(get_star_alpha());

	if (!have_sun) {
		bkg_color = BACKGROUND_NIGHT;
	}
	else {
		blend_color(bkg_color, BACKGROUND_NIGHT, attenuate_sun_color(base_sky_color), star_alpha, 1);
	}
	if (is_cloudy) {
		colorRGBA const orig_bkgc(bkg_color);
		blend_color(bkg_color, bkg_color, GRAY, 0.5, 1);
		UNROLL_3X(bkg_color[i_] = min(bkg_color[i_], orig_bkgc[i_]);) // can't make it brighter
	}
	//if (atmosphere < 1.0) {blend_color(bkg_color, bkg_color, BACKGROUND_NIGHT, atmosphere, 0);}
}


void add_sun_effect(colorRGBA &color) {

	color = color.modulate_with(sun_color);
	float cmult(have_sun ? get_lf_scale(light_factor) : 0.0); // light from sun

	if (!combined_gu && light_factor < 0.6) { // light from moon
		float const lfs(get_lf_scale(get_moon_light_factor()));
		cmult += 0.5*CLIP_TO_01(5.0f*(0.6f - light_factor))*lfs;
	}
	color *= (0.8*cmult + 0.2);
}


void setup_linear_fog(colorRGBA const &color, float fog_end) {

	cur_fog_color = color;
	cur_fog_end   = fog_end;
	add_sun_effect(cur_fog_color);
}


void auto_advance_camera() {

	if (run_forward && !advanced) advance_camera(MOVE_FRONT);
	advanced = 0;
}


void config_bkg_color_and_clear(bool no_fog) {

	calc_bkg_color();
	glClearColor_rgba((!no_fog && show_fog) ? GRAY : bkg_color);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); // Clear the background
}


void reset_planet_defaults() {

	have_sun   = 1;
	atmosphere = def_atmosphere;
	vegetation = def_vegetation;
	base_cloud_color = WHITE;
	base_sky_color   = BACKGROUND_DAY;
}


float get_light_pos_scale() {
	return ((world_mode == WMODE_INF_TERRAIN) ? 10.0 : 1.0); // hack: make the sun and moon far away in inf terrain mode 
}

void setup_sun_moon_light_pos() {
	set_gl_light_pos(0, sun_pos *get_light_pos_scale(), LIGHT_W_VAL);
	set_gl_light_pos(1, moon_pos*get_light_pos_scale(), LIGHT_W_VAL);
}

void setup_lighting() {
	
	// background color code
	config_bkg_color_and_clear(world_mode == WMODE_INF_TERRAIN);

	// lighting code - RGB intensity for ambient and diffuse (specular is set elsewhere per object)
	float const mlf(get_moon_light_factor());
	colorRGBA ambient, diffuse;
	ambient[3] = diffuse[3] = 1.0;
	enabled_lights = 0;

	// Note: should this be set in universe lighting?
	colorRGBA ambient_c;
	blend_color(ambient_c, bkg_color, attenuate_sun_color(WHITE), 0.5, 1);

	for (unsigned i = 0; i < 3; ++i) {
		diffuse[i] = DSCALE*sunlight_color[i];
		ambient[i] = ASCALE*ambient_c[i];
	}
	for (unsigned i = 0; i < 3; ++i) {
		if (is_cloudy) {
			diffuse[i] -= (auto_time_adv ? 0.06  : 0.15);
			ambient[i] -= (auto_time_adv ? 0.025 : 0.06);
		}
		diffuse[i] -= 0.12*cloud_cover;
		ambient[i] -= 0.05*cloud_cover;
		ambient[i] *= ambient_scale;
	}
	if (light_factor <= 0.4) { // moon
		calc_moon_atten(ambient, diffuse, mlf);
		set_colors_and_enable_light(1, ambient, diffuse*moonlight_brightness);
	}
	else if (light_factor >= 0.6) { // sun
		for (unsigned i = 0; i < 3; ++i) { // add more brightness
			diffuse[i] += 0.2;
			ambient[i] += 0.1;
		}
		set_colors_and_enable_light(0, ambient, diffuse*sunlight_brightness);
	}
	else { // sun and moon
		float const lfd(get_light_factor_scale()), lfn(1.0 - lfd);

		for (unsigned i = 0; i < 3; ++i) { // should diffuse depend more on angle than ambient?
			diffuse[i] = (diffuse[i] + 0.2)*lfd;
			ambient[i] = (ambient[i] + 0.1)*lfd;
		}
		set_colors_and_enable_light(0, ambient, diffuse*sunlight_brightness); // sun

		for (unsigned i = 0; i < 3; ++i) {
			diffuse[i] = (diffuse[i]/lfd - 0.2)*lfn;
			ambient[i] = (ambient[i]/lfd - 0.1)*lfn;
		}
		calc_moon_atten(ambient, diffuse, mlf);
		set_colors_and_enable_light(1, ambient, diffuse*moonlight_brightness); // moon
	}

	// setup light position (after enabling lights)
	setup_sun_moon_light_pos();
	setup_gl_light_atten(0, 1.0, 0.0, 0.0); // reset attenuation to 1.0
}


void draw_sun_moon_stars(bool no_update) {

	float star_alpha(get_star_alpha(is_cloudy != 0));
	if (star_alpha > 0.0) {gen_and_draw_stars(star_alpha, 0, no_update);}
	if (light_factor <= 0.6 || atmosphere <= 0.5) {draw_moon();} // moon
	if (light_factor >= 0.4) {draw_sun();} // sun
}


void draw_universe_bkg() {

	RESET_TIME;
	// setup sun position and related parameters
	update_sun_and_moon();
	point const camera(get_camera_pos()), new_sp(sun_pos.get_norm()), player_pos(get_player_pos());
	init_universe_display();
	sun_pos = univ_sun_pos - player_pos; // univ_sun_pos won't be correct on the first call
	point const old_sp(sun_pos.get_norm());
	rotate_vector3d_by_vr(old_sp, new_sp, sun_pos);
	sun_pos /= UNIV_NCLIP_SCALE;
	set_player_pdu(new_sp, old_sp);
	do_look_at();

	// transform universe coordinate system into current mesh coordinate system
	fgPushMatrix();
	point camera_r(camera);
	rotate_from_v2v(new_sp, old_sp);
	rotate_vector3d_by_vr(new_sp, old_sp, camera_r);
	translate_to(camera_r - player_pos); // shift player_pos to camera pos

	// draw universe as background
	config_bkg_color_and_clear(1);
	point const camera_pos_orig(camera_pos);
	camera_pos = player_pos; // trick universe code into thinking the camera is at the player's ship
	stop_player_ship();
	if (TIMETEST) PRINT_TIME("0.1");
	bool const camera_above_clouds(world_mode == WMODE_INF_TERRAIN && camera_pos_orig.z > get_tt_cloud_level());
	bool const no_stars((is_cloudy || (atmosphere > 0.8 && light_factor >= 0.6)) && !camera_above_clouds), no_asteroid_dust(no_stars);
	draw_universe(1, 1, 1, (no_stars ? 2 : 0), 0, no_asteroid_dust); // could clip by horizon?
	if (TIMETEST) PRINT_TIME("0.2");
	camera_pos = camera_pos_orig;
	fgPopMatrix(); // undo universe transform

	// setup sun light source
	float const sun_intensity(max(0.25f, min(4.0f, 1000.0f*univ_sun_rad/sun_pos.mag())));
	setup_current_system(sun_intensity);
	set_gl_light_pos(0, sun_pos*get_light_pos_scale(), LIGHT_W_VAL);
	disable_light(1); // no moonlight (for now)
	if (!have_sun || light_factor < 0.5) {set_light_ds_color(0, BLACK);} // sun below horizon: no diffuse or specular
	check_zoom(); // reset perspective
	//light_factor = (PI + 2.0*asinf(dot_product(sun_pos.get_norm(), plus_z)))/TWO_PI; // shouldn't change much from previous light_factor

	// setup background and init for standard mesh draw
	if (light_factor > 0.4) { // translucent blue for atmosphere
		colorRGBA color(bkg_color);
		color.alpha *= 0.75*atmosphere*min(1.0, (light_factor - 0.4)/0.2);
		vector<camera_filter> cfs;
		cfs.push_back(camera_filter(color, 1, -1, 0));
		draw_camera_filters(cfs);
	}
	do_look_at();
	glClear(GL_DEPTH_BUFFER_BIT);
	if (TIMETEST) PRINT_TIME("0.3");
}


void draw_game_elements(int timer1, int reflection_pass=0) {

	if (TIMETEST) PRINT_TIME("U");
	draw_camera_weapon(1, reflection_pass);
	draw_projectile_effects(reflection_pass);
	if (TIMETEST) PRINT_TIME("V");
	draw_smoke_and_fires();
	if (TIMETEST) PRINT_TIME("W");
	draw_scheduled_weapons(reflection_pass == 0);
}


void setup_basic_fog() {
	if (show_fog) {reset_fog();	fog_enabled = 1;}
}


void add_uw_light_color_comp(point const &lpos, float weight, colorRGBA &color) {

	// check for in shadow? what about tiled terrain?
	weight *= 0.5 + 0.5*max(0.0f, lpos.z/lpos.mag()); // vertical component (which penetrates water)
	UNROLL_3X(color[i_] += weight;)
}


void atten_uw_fog_color(colorRGBA &color, float depth) {

	colorRGBA light_color(BLACK);
	float const lf(get_lf_scale(light_factor));
	if (lf > 0.0) {add_uw_light_color_comp(sun_pos,  lf,           light_color);}
	if (lf < 1.0) {add_uw_light_color_comp(moon_pos, 0.5*(1.0-lf), light_color);}
	if (is_cloudy) {light_color *= 0.5;}
	color  = color.modulate_with(light_color);
	atten_by_water_depth(&color.R, depth);
	color *= FOG_COLOR_ATTEN;
	colorRGBA filt_color(color);
	filt_color.A = 0.25;
	add_camera_filter(filt_color, 1, -1, CAM_FILT_UWATER);
}


void set_inf_terrain_fog(bool underwater, float zmin2) {

	float fog_dist;
	colorRGBA fog_color;

	if (underwater) { // under water/ice
		float const camera_z(get_camera_pos().z);
		fog_color = colorRGBA(get_tt_water_color(), 1.0); // alpha = 1.0
		atten_uw_fog_color(fog_color, 2.0f*water_params.alpha*(water_plane_z - camera_z)); // more opaque = effectively deeper
		fog_dist = (water_is_lava ? 0.5f : 1.0f)*(0.3f + 1.5f*Z_SCENE_SIZE*(camera_z - zmin2)/max(1.0E-3f, (water_plane_z - zmin2))) * max(0.1f, (1.5f - water_params.alpha));
	}
	else {
		get_avg_sky_color(fog_color);
		colorRGBA const cloud_color(get_cloud_color(), 1.0); // alpha = 1.0
		blend_color(fog_color, cloud_color, bkg_color, 0.375, 1); // weighted more towards bkg_color
		if (water_is_lava) {blend_color(fog_color, fog_color, colorRGBA(1.0, 0.2, 0.0, 1.0), 0.75, 0);} // add slight red-orange color
		fog_dist = get_inf_terrain_fog_dist();
		if      (player_in_tunnel  ) {fog_color = BLACK     ;} // black fog in tunnels
		else if (player_in_basement) {fog_color = GRAY_BLACK;} // dark gray fog in basements
	}
	setup_linear_fog(fog_color, fog_dist); // under water/ice
	fog_enabled = 1;
}


void scroll_scene() {

	RESET_TIME;
	point const camera(get_camera_pos());
	cout << endl << "Shifting + "; // will produce "Shifting + Generating scene"
	camera_change = 1;
	scrolling     = 1;
	dx_scroll     = int(camera.x*DX_VAL_INV);
	dy_scroll     = int(camera.y*DY_VAL_INV);
	vector3d const vd(-DX_VAL*dx_scroll, -DY_VAL*dy_scroll, 0.0);
	surface_pos  += vd;
	xoff2        += dx_scroll;
	yoff2        += dy_scroll;
	shift_all_objs(vd);
	gen_scene(1, 1, 1, 0, 0);
	regen_lightmap(); // not shiftable
	if (display_mode & 0x04) {water_plane_z = get_water_z_height();}
	update_temperature(0);
	invalidate_snow_coverage();
	recreated = 1;
	scrolling = 0;
	clear_landscape_vbo = 1;
	PRINT_TIME("*** Top Level: Final");
}


float get_ocean_wave_height() {

	if (water_is_lava || !(display_mode & 0x0100)) return 0.0;
	static float time(0.0);
	if (animate2 && temperature > W_FREEZE_POINT) {time += fticks;}
	return ocean_wave_height*sin(1.0*time/TICKS_PER_SECOND); // add small waves
}


void draw_sun_flare(int ignore_cobj=-1, float intensity=1.0) {

	//RESET_TIME;
	if (!is_sun_flare_visible()) return;
	point const sun_pos(get_sun_pos()), viewer(get_camera_pos());
	//if (is_cloudy) {intensity *= 0.7;}

	if (world_mode == WMODE_GROUND) {
		unsigned const npts = 16;
		static point pts[npts];
		static bool pts_valid(0);
		float tot_light(0.0);
		vector3d const view_dir((sun_pos - viewer).get_norm());
		point viewer_pos(viewer);
		if (enable_clip_plane_z && view_dir.z > 0.0 && viewer.z < clip_plane_z) {viewer_pos += view_dir*((clip_plane_z - viewer.z)/view_dir.z);} // move above clip plane
		
		for (unsigned i = 0; i < npts; ++i) {
			int index; // unused
			if (!pts_valid) {pts[i] = signed_rand_vector_norm();}
			point const pos(sun_pos + pts[i]*sun_radius);

			if (coll_pt_vis_test(pos, viewer_pos, 0.0, index, ignore_cobj, 0, 1) && (!(display_mode & 0x01) || !line_intersect_mesh(pos, viewer_pos, 0))) {
				tot_light += 1.0 - get_cloud_density(viewer_pos, view_dir);
			}
		}
		pts_valid = 1;
		if (tot_light == 0) return;
		intensity *= 0.1 + 0.9*tot_light/npts;
		if (show_fog)  {intensity *= 0.4;}
	}
	else if (world_mode == WMODE_INF_TERRAIN) {
		if (sun_pos.z < zmin) return; // sun below the mesh
		if (viewer.z < water_plane_z) {intensity *= CLIP_TO_01(1.0f - 1.0f*(water_plane_z - viewer.z));} // attenuate sun flare when underwater
	}
	DoFlares(viewer, camera_origin, sun_pos, 1.0, (combined_gu ? 15.0*univ_sun_rad : 1.0), intensity);
	//PRINT_TIME("Query + Flare");
}


void set_multisample(bool enable) {
	if (!enable_multisample) return;
	if (enable) {glEnable(GL_MULTISAMPLE);} else {glDisable(GL_MULTISAMPLE);}
}

void setup_depth_clamp() {
	if (enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);} else {glDisable(GL_DEPTH_CLAMP);}
}

void draw_sky_and_clouds(bool camera_side, bool no_update=0) {
	draw_sky(camera_side, no_update);
	draw_puffy_clouds(camera_side, no_update);
}

void create_reflection_and_portal_textures() {
	//timer_t timer("Create Reflection Textures"); // 340 => 250 => 230 => 195 => 184 => 112 => 103 => 98 => 92 => 82 => 72 => 63 => 54 => 47
	if (enable_reflection_plane()) {create_gm_z_reflection();} // must be before draw background but after setup_object_render_data()
	if (!enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);} // enable depth clamp if not yet enabled - useful for cube maps
	ensure_model_reflection_cube_maps();
	reflective_cobjs.create_textures();
	if (!enable_depth_clamp) {glDisable(GL_DEPTH_CLAMP);} // restore orig value if needed
	create_portal_textures();
}

void flashlight_next_frame() {
	static bool last_flashlight_on(0);
	if (flashlight_on != last_flashlight_on) {play_switch_weapon_sound();}
	last_flashlight_on = flashlight_on;
	flashlight_on = 0;
}

void maybe_update_loading_screen(const char *str) {
	if (!in_loading_screen) return;
	if (shader_is_active()) return; // don't update loading screen in the middle of drawing something else
	static string suffix;
	check_gl_error(566);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	fgPushMatrix();
	fgLoadIdentity();
	draw_text(PURPLE, -0.01, 0.0, -0.02, (string("Loading: ") + str + suffix), 1.0);
	fgPopMatrix();
	glutSwapBuffers();
	bool const had_error(check_gl_error(567));
	if (had_error) {in_loading_screen = 0;}
	suffix.push_back('.');
}

void begin_loading_screen() {
	static bool showed_loading_screen(0);
	if (showed_loading_screen) return; // already showed screen
	showed_loading_screen = 1;
	config_bkg_color_and_clear(1);
	bind_vao(0); // set to default VAO
	draw_text(PURPLE, 0.0, 0.008, -0.02, "Loading...", 2.5);
	swap_buffers_and_redraw();
	in_loading_screen = (world_mode == WMODE_INF_TERRAIN); // only in tiled terrain mode for now (for slow city/buildings scene)
}

void display() {

	check_gl_error(0);
	static unsigned counter(0);
	//if (counter <= 120) {cout << "frame " << counter << " time " << GET_TIME_MS() << "ms" << endl;} // TESTING: 26s for config_heightmap with people+cars
	//print_texture_memory_usage(); // TESTING

	// hack to avoid slow frames when starting OpenMP threads; for some reason, threads started during the first 100 frames slow those frames down
	if (counter++ < 100) {
		glutSwapBuffers();
		post_window_redisplay();
		check_gl_error(11);
		return;
	}
	if (start_maximized) {
		toggle_fullscreen();
		begin_loading_screen();
		start_maximized = 0;
		return;
	}
	RESET_TIME;
	static int init(0), tticks(0);
	static point old_spos(0.0, 0.0, 0.0);
	++cur_display_iter;
	num_frame_draw_calls = 0;
	proc_kbd_events();

	if (!init) { // the first frame
		init   = 1;
		fticks = 1.0;
		time0  = timer1;
		begin_loading_screen(); // for the !start_maximized case
	}
	else if (animate && !DETERMINISTIC_TIME) {
		double ftick(0.0);
		static float carry(0.0);
		double const time_delta((TICKS_PER_SECOND*(timer1 - time0))/1000.0f);

		if (reset_timing) {
			iticks  = 0;
			ftick   = 0.0;
			carry   = 0.0;
			tticks  = int(time_delta);
		}
		else {
			assert(timer1 >= time0);
			tfticks = time_delta;
			ftick   = (time_delta - tticks);
			iticks  = (int)ftick;
			tticks += iticks;
			carry   = ftick - (double)iticks;
			ftick   = min(ftick, 20.0);
			if (animate2) {sim_ticks = tfticks;}
		}
		fticks = max(TOLERANCE, float(0.9*fticks + 0.1*(ftick - carry))); // slow averaging filter
		assert(fticks >  0.0 && fticks < 1.0E12);
		assert(iticks >= 0   && iticks < 1000000000);
	}
	else {
		fticks = 1.0;
		iticks = 1;
	}
	flashlight_next_frame();
	tstep         = TIMESTEP*fticks;
	reset_timing  = 0;
	check_gl_error(1);
	set_fill_mode();
	set_std_blend_mode();

	if (map_mode && world_mode != WMODE_UNIVERSE) {
		if (!pause_frame) {uevent_advance_frame();}

		if (world_mode == WMODE_INF_TERRAIN) { // map mode infinite terrain
			camera_origin = surface_pos;
			update_cpos();
			apply_camera_offsets(get_camera_pos());
			check_xy_offsets();
			next_city_frame(0); // make sure the cars animate
		}
		else if (world_mode == WMODE_GROUND) {
			process_groups();
			build_cobj_tree(1, 0); // ensure smiley positions are valid for map mode (since framrate is low)
			if (game_mode) {update_blasts(); update_game_frame();}
		}
		draw_overhead_map();
		if (world_mode == WMODE_INF_TERRAIN) {show_crosshair(WHITE, 0);}
		swap_buffers_and_redraw();
		return;
	}
	up_vector.assign(0.0, sinf(up_theta), camera_y*cosf(up_theta));
	setup_sphere_vbos();
	check_gl_error(2);
	update_sound_loops();
	set_multisample(1);
	glEnable(GL_DEPTH_TEST); // Note: seems to be required to make reflections work after first un-maximize (why?)
	setup_depth_clamp();
	bind_vao(0); // set to default VAO
#ifdef USE_GPU_TIMER
	gpu_timer_t gpu_timer;
#endif
	
	if (world_mode == WMODE_UNIVERSE) {
		in_loading_screen = 0; // if we got here, loading is done
		display_universe(); // infinite universe
	}
	else {
		if (!pause_frame) {uevent_advance_frame();}
		earth_radius = 2.0;
		sun_radius   = 1.5;
		moon_radius  = 2.0;
		check_zoom();

		if (init_x && world_mode == WMODE_GROUND) {
			gen_scene(1, 1, KEEP_MESH, 0, 0);
			init_x   = 0;
			show_fog = 0;
		}
		
		// timing and framerate code
		float const framerate(framerate_tracker.show_cur_framerate());
		if (world_mode == WMODE_GROUND) {process_platforms_falling_moving_and_light_triggers();} // must be before camera code
		else if (world_mode == WMODE_INF_TERRAIN) {camera_mode = 1;} // force to ground/walking mode
		if (world_mode == WMODE_GROUND) {UNLIMITED_WEAPONS = config_unlimited_weapons;}
		else if (world_mode == WMODE_INF_TERRAIN) {UNLIMITED_WEAPONS = 1;} // always have unlimited weapons in tiled terrain mode

		// camera position code
		auto_advance_camera();
		if (world_mode == WMODE_GROUND && camera_mode == 1 && !spectate) {check_popup_text();}

		if (camera_mode == 1 && camera_surf_collide && !spectate) {
			force_onto_surface_mesh(surface_pos);

			if (world_mode == WMODE_INF_TERRAIN && temperature <= W_FREEZE_POINT) {
				surface_pos.z = max(surface_pos.z, (water_plane_z + CAMERA_RADIUS)); // camera on ice in WM3
			}
		}
		else {
			remove_reset_coll_obj(camera_coll_id);
		}
		update_camera_velocity(surface_pos - old_spos);
		old_spos = surface_pos;
		static double temp_c_radius(0.0);

		if (camera_mode == 1 && !camera_view) {
			cpos2       = surface_pos;
			camera_view = 1;
		}
		if (camera_view) {
			if (c_radius >= C_RADIUS0) {temp_c_radius = c_radius;}
			c_radius      = (double)CR_SCALE*C_RADIUS0;
			camera_origin = cpos2;
			if (camera_mode == 1 && !spectate) {camera_origin.z += get_player_height();}
		}
		else if (temp_c_radius >= C_RADIUS0) {
			c_radius      = temp_c_radius;
			temp_c_radius = 0.0;
		}
		update_cpos();
		point const camera(get_camera_pos());
		float depth(0.0); // used for underwater fog calculation
		underwater = check_underwater(CAMERA_ID, depth);
		auto_advance_time();
		if (animate2) {total_wind += wind*fticks;}
		check_gl_error(3);
		if (TIMETEST) PRINT_TIME("\n\nA");
		
		if (!combined_gu) {
			do_look_at();
			sun_color = sunlight_color;
			apply_red_sky(sun_color);
			reset_planet_defaults();
			setup_lighting();
			check_gl_error(4);
			if (TIMETEST) PRINT_TIME("B");
		}
		if (world_mode == WMODE_INF_TERRAIN) { // infinite terrain mode
			display_inf_terrain();
		}
		else { // finite terrain mode
			if (mesh_invalidated) {
				gen_mesh_bsp_tree();
				clear_landscape_vbo = 1;
				mesh_invalidated    = 0;
			}
			player_in_water = 0; // only used for TT mode
			if (show_fog || underwater) {fog_enabled = 1;}
			if (!show_lightning && animate2) {l_strike.disable();} // Note: only legal to call this when lighting is updated (animate2)
			compute_brightness();
			if (TIMETEST) PRINT_TIME("D");

			// run physics and collision detection
			process_groups();
			check_gl_error(12);
			if (TIMETEST) PRINT_TIME("E");
			if (b2down) {fire_weapon();}
			update_weapon_cobjs(); // and update cblade
			setup_dynamic_teleporters();
			check_gl_error(6);
			proc_voxel_updates(); // with the update here, we avoid making the voxels and shadows out of sync
			if (TIMETEST) PRINT_TIME("F");

			//proc_voxel_updates(); // with the update here, we avoid uploading the modified voxel VBOs during shadow map rendering

			// send data to GPU
			setup_object_render_data();
			check_gl_error(101);
			in_loading_screen = 0; // if we got here, loading is done

			// create shadow map
			if (combined_gu) {do_look_at();}
			create_shadow_map(); // where should this go? must be after draw_universe_bkg()
			if (TIMETEST) PRINT_TIME("G");

			create_reflection_and_portal_textures();

			// draw background
			if (combined_gu) {draw_universe_bkg();} // infinite universe as background
			else {
				draw_sun_moon_stars(0);
				draw_earth();
			}
			draw_sky_and_clouds(0); // Note: depth test is disabled, so must be drawn first
			check_gl_error(5);
			if (TIMETEST) PRINT_TIME("G2");

			if (underwater) {
				colorRGBA fog_color(((water_is_lava ? LAVA_COLOR : (temperature <= W_FREEZE_POINT) ? ICE_C : WATER_C)), 1.0); // under ice/water, alpha = 1.0
				select_liquid_color(fog_color, camera);
				atten_uw_fog_color(fog_color, depth);
				float const fog_dist(0.2f + (0.25f + 0.75f*fog_color.B)*(1.5f*Z_SCENE_SIZE)*(camera.z - zmin)/((camera.z + depth) - zmin));
				setup_linear_fog(fog_color, (water_is_lava ? 0.5 : 1.0)*fog_dist);
			}

			// draw the scene
			draw_camera_weapon(0);
			if (TIMETEST) PRINT_TIME("H");

			draw_coll_surfaces(0, 0);
			if (TIMETEST) PRINT_TIME("I");
			
			if (display_mode & 0x01) {display_mesh();} // draw mesh
			check_gl_error(7);
			if (TIMETEST) PRINT_TIME("J");

			draw_grass();
			draw_scenery();
			if (TIMETEST) PRINT_TIME("K");

			draw_solid_object_groups();
			check_gl_error(8);
			if (TIMETEST) PRINT_TIME("L");

			draw_stuff(!underwater, timer1);
			if (TIMETEST) PRINT_TIME("M");

			if (display_mode & 0x04) {draw_water();} // must be after process_groups()
			check_gl_error(9);
			if (TIMETEST) PRINT_TIME("N");
			draw_stuff(underwater, timer1);
			if (TIMETEST) PRINT_TIME("T");

			if (show_lightning && l_strike.is_enabled() && animate2) {
				if ((rand()&1) == 0) {gen_smoke(l_strike.get_hit_pos());}
				if ((rand()&7) == 0) {gen_fire(l_strike.get_hit_pos(), 1.0, NO_SOURCE);}
			}
			update_blasts(); // not really an update, but needed for draw_blasts
			draw_game_elements(timer1);
			setup_basic_fog();
			draw_sky_and_clouds(1);
			draw_sun_flare(camera_coll_id);
			check_gl_error(10);
			//draw_scene_bounds_and_light_frustum(get_light_pos()); // TESTING
		} // WMODE_GROUND
		check_gl_error(13);
		if (debug_event_pos != all_zeros) {show_debug_event_pos(debug_event_pos);}
		final_draw(framerate);
		purge_coll_freed(0); // optional
		camera_flight = 0;
		if (game_mode) {update_game_frame();} // even in TT mode

		if (game_mode && world_mode == WMODE_GROUND) {
			show_user_stats();
			show_blood_on_camera();
			show_crosshair(WHITE, do_zoom);
			inf_terrain_fire_mode = 0;
		}
		else if (world_mode == WMODE_INF_TERRAIN && (inf_terrain_fire_mode || tt_fire_button_down)) {
			show_crosshair(get_inf_terrain_mod_color(), do_zoom);
		}
		else if (world_mode == WMODE_INF_TERRAIN && show_bldg_pickup_crosshair) {
			show_crosshair((can_pickup_bldg_obj ? ((can_pickup_bldg_obj == 2) ? RED : GREEN) : (can_do_building_action ? CYAN : WHITE)), can_pickup_bldg_obj);
		}
		else if (world_mode == WMODE_GROUND && voxel_editing) {
			show_crosshair(((voxel_editing == 2) ? RED : GREEN), do_zoom);
		}
		else if (spraypaint_mode) {
			draw_spraypaint_crosshair();
		}
		if (display_framerate && !game_mode && ((world_mode == WMODE_INF_TERRAIN && !show_bldg_pickup_crosshair) || (world_mode == WMODE_GROUND && camera_mode == 1))) {
			draw_compass_and_alt();
		}
		if (world_mode == WMODE_GROUND && camera_mode == 1 && camera_surf_collide) {show_player_keycards();}
		if (world_mode == WMODE_INF_TERRAIN) {building_gameplay_next_frame();} // Note: must be after crosshair drawing
		
		if (display_framerate) { // notify the user of lighting progress
			if (indir_lighting_updated()) {draw_text(PURPLE, 0.007*(float)window_width/(float)window_height, -0.009, -0.02, "Lighting Updating");}
			else if (!lighting_update_text.empty()) {
				draw_text(PURPLE, 0.007*(float)window_width/(float)window_height, -0.009, -0.02, lighting_update_text);
				lighting_update_text.clear();
			}
		}
		tt_fire_button_down = 0;
		if (TIMETEST) PRINT_TIME("X");

		if (dynamic_mesh_scroll && world_mode == WMODE_GROUND && camera_mode == 1 && !camera_view) {
			float const cdist(max(fabs(camera.x/X_SCENE_SIZE), fabs(camera.y/Y_SCENE_SIZE)));
			if (cdist > REL_SCROLL_DIST) {scroll_scene();}
		}
	} // not universe mode
#ifdef USE_GPU_TIMER
	gpu_timer.show();
#endif
	draw_enabled_ui_menus();
	swap_buffers_and_redraw();
	check_gl_error(11);
	if (TIMETEST) PRINT_TIME("Y");
}


void display_universe() { // infinite universe

	static int init(0);
	RESET_TIME;

	if (!init || init_x) {
		init     = 1;
		init_x   = 0;
		show_fog = 0;
		update_cpos();
	}
	camera_view = 0;
	float const framerate(framerate_tracker.get_framerate());
	setup_landscape_tex_colors(ALPHA0, ALPHA0); // reset for asteroids
	init_universe_display();
	bkg_color   = BACKGROUND_NIGHT;
	glClearColor_rgba(bkg_color);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); // Clear the background

	check_gl_error(30);
	auto_advance_camera();
	if (TIMETEST) PRINT_TIME("\nSetup");
	apply_univ_physics(); // physics loop
	if (TIMETEST) PRINT_TIME("Physics");
	int const last_csc(camera_surf_collide);
	camera_mode         = 1;
	camera_surf_collide = 0;
	c_radius            = C_RADIUS0;
	camera_origin       = get_player_pos();
	update_cpos();
	if (!pause_frame) {uevent_advance_frame();}
	do_look_at();
	if (b2down) {fire_weapon();} // just sets fire_key=1
	check_gl_error(31);
	draw_universe();
	check_gl_error(32);
	if (TIMETEST) PRINT_TIME("Draw Universe");
	draw_universe_blasts();
	if (TIMETEST) PRINT_TIME("Draw Blasts");
	final_draw(framerate);
	show_crosshair(WHITE, do_zoom);
	if (display_framerate) {draw_universe_stats();}
	camera_surf_collide = last_csc;
	check_gl_error(33);
}


// Note: assumes the camera is not underwater
void draw_scene_from_custom_frustum(pos_dir_up const &pdu, int cobj_id, int reflection_pass, bool inc_mesh, bool inc_grass, bool inc_water) {

	check_gl_error(532);
	pos_dir_up const prev_camera_pdu(camera_pdu);
	point const prev_camera_pos(camera_pos);
	camera_pdu = pdu;
	camera_pos = pdu.pos;
	bool const disable_occ_cull(reflection_pass == 1 && (display_mode & 0x08) != 0); // disable occlusion culling
	if (disable_occ_cull) {display_mode &= ~0x08;}

	// draw background
	if (!combined_gu) { // Note: universe not drawn here (must be drawn by caller before setting up matrices)
		draw_sun_moon_stars(1);
		draw_earth();
	}
	draw_sky_and_clouds(0, 1); // Note: depth test is disabled, so must be drawn first
	
	// draw the scene
	draw_camera_weapon(0, reflection_pass);
	update_shadow_matrices();
	draw_coll_surfaces(0, reflection_pass);
	
	// the mesh and grass are generally under the reflection plane, so can be skipped in the reflection plane pass
	if (inc_mesh && (display_mode & 0x01)) {display_mesh(0, 1);} // draw mesh
	if (inc_grass) {draw_grass();}
	draw_scenery();
	draw_solid_object_groups(reflection_pass);
	draw_stuff(1, 0, reflection_pass);
	if (inc_water && (display_mode & 0x04)) {draw_water(1, (reflection_pass == 2));} // if we can exclude the mesh we can probably exclude the water as well
	draw_stuff(0, 0, reflection_pass);
	draw_game_elements(0, reflection_pass);
	setup_basic_fog();
	draw_sky_and_clouds(1);
	draw_sun_flare(cobj_id);

	if (camera_mode == 1 && camera_surf_collide && begin_motion) { // player is on the ground and collision in enabled, draw the player smiley
		draw_player_model(surface_pos + vector3d(0.0, 0.0, get_player_height()), calc_camera_direction(), int(tfticks));
	}
	// restore original values
	if (disable_occ_cull) {display_mode |= 0x08;}
	camera_pdu = prev_camera_pdu;
	camera_pos = prev_camera_pos;
}


void run_tt_gameplay() {

	process_groups();
	setup_dynamic_teleporters();
	draw_select_groups(1);
	draw_select_groups(0);
	update_blasts();
}


void draw_tiled_terrain_and_transparent_geom(float terrain_zmin, unsigned tt_reflection_tid, bool draw_water, bool camera_above_clouds) {
	draw_tiled_terrain(0);
	render_tt_models(0, 1); // transparent pass
	//if (underwater ) {draw_local_precipitation();}
	if (draw_water ) {draw_water_plane(water_plane_z, terrain_zmin, tt_reflection_tid);}
	if (show_lightning) {end_tiled_terrain_lightning();}

	if (!underwater) {
		maybe_draw_rainbow();
		draw_tiled_terrain_clouds(0);
		draw_local_precipitation();
	}
	else {
		draw_underwater_particles(terrain_zmin);
	}
	draw_cloud_planes(terrain_zmin, 0, camera_above_clouds, 0);
}

void display_inf_terrain() { // infinite terrain mode (Note: uses light params from ground mode)

	static int init_xx(1);
	RESET_TIME;
	//timer_t timer("Display Inf Terrain"); // 6.9 no update / 10.6 1-thread / 8.0 2-threads / 7.6 3-threads

	if (init_x || init_xx) {
		init_xx  = 0;
		show_fog = 1;
		c_radius = 0.001;
		update_cpos();
		if (!player_custom_start_pos) {surface_pos = all_zeros;}
		camera_surf_collide = 1;
	}
	camera_view = 0;
	if (camera_surf_collide) {check_player_tiled_terrain_collision();}
	follow_city_actor(); // after collision detection so that it doesn't apply to the actor we're following
	update_temperature(0);
	apply_camera_offsets(get_camera_pos());
	compute_brightness();
	set_global_state();
	if (b2down) {fire_weapon();}
	update_weapon_cobjs(); // and update cblade

	// drawing
	bool const water_enabled((display_mode & 0x04) && !DISABLE_WATER);
	water_plane_z = (water_enabled ? get_max_sea_level() : -10*FAR_DISTANCE);
	camera_mode   = 1; // walking on ground
	float min_camera_dist(0.0);
	float const terrain_zmin(update_tiled_terrain(min_camera_dist));
	bool const change_near_far_clip(!camera_surf_collide && min_camera_dist > 0.0 && !do_zoom);
	bool const draw_water(water_enabled && water_plane_z >= terrain_zmin && !player_in_basement); // not correct if the basement is below the water line
	if (show_fog || underwater) {set_inf_terrain_fog(underwater, terrain_zmin);}
	unsigned tt_reflection_tid(0);
	if (draw_water && !underwater) {tt_reflection_tid = create_tt_reflection(terrain_zmin);}
	create_building_reflections();
	far_clip_ratio = 1.0; // reset to default, may be overwritten below

	if (enable_tt_model_indir) {
		calc_cur_ambient_diffuse(); // required to handle lighting updates
		upload_smoke_indir_texture();
	}
	if (combined_gu) {
		draw_universe_bkg(); // infinite universe as background
		check_gl_error(4);
	}
	else {
		config_bkg_color_and_clear(1);
		draw_sun_moon_stars(0);
	}
	if (change_near_far_clip) {
		float const near_clip(NEAR_CLIP + 0.01*min_camera_dist);
		float const far_clip(get_tt_fog_based_far_clip(min_camera_dist));
		set_perspective_near_far(near_clip, far_clip);
		do_look_at(); // clear depth buffer?
		camera_pdu.near_ = near_clip; // override camera frustum near/far clip so that VFC will be correct
		camera_pdu.far_  = far_clip;
		far_clip_ratio   = far_clip/FAR_CLIP;
	}
	//draw_puffy_clouds(0);
	draw_camera_weapon(0);
	bool const camera_above_clouds(get_camera_pos().z > get_tt_cloud_level());
	draw_cloud_planes(terrain_zmin, 0, !camera_above_clouds, 1); // these two lines could go in either order
	draw_sun_flare();
	if (TIMETEST) PRINT_TIME("3.2");
	if (show_lightning) {draw_tiled_terrain_lightning(0);}
	// load all tiles/buildings around the player; since there's a cap on data generated per frame, we use the first 20 frames
	bool const disable_vfc(pre_load_full_tiled_terrain && frame_counter <= 20), prev_pdu_valid(camera_pdu.valid);
	if (disable_vfc) {camera_pdu.valid = 0;} // disable view frustum culling
	pre_draw_tiled_terrain(); // must be before render_tt_models()
	if (in_loading_screen) {glutSwapBuffers();} // show our cloud background
	in_loading_screen = 0; // if we got here, loading is done
	if (TIMETEST) PRINT_TIME("3.26");
	render_tt_models(0, 0); // opaque pass; draws city buildings, cars, etc.
	if (TIMETEST) PRINT_TIME("3.27");

	// threads: 0=draw tiled terrain (2.5ms) + transparent (0.3), 1=update roads and cars (2.3ms), 2=pedestrians (3.8ms) + building AI (0.85ms)
	// Note: it's questionable to update (move) cars between the opaque and transparent pass because the parts will be out of sync;
	// however, only the headlight flares are drawn in the transparent pass, and it doesn't seem to be a problem, so we allow it
	if (have_city_models()) {
		//timer_t timer("City Update MT"); // 4.65ms
#pragma omp parallel num_threads(3)
		if (omp_get_thread_num_3dw() == 0) {draw_tiled_terrain_and_transparent_geom(terrain_zmin, tt_reflection_tid, draw_water, camera_above_clouds);} // drawing must be on thread 0
		else {next_city_frame(1);} // other threads (if threads enabled, else serial)
	}
	else { // serial version
		//timer_t timer("City Update"); // 10.0ms
		next_city_frame(0);
		draw_tiled_terrain_and_transparent_geom(terrain_zmin, tt_reflection_tid, draw_water, camera_above_clouds);
	}
	camera_pdu.valid = prev_pdu_valid; // restore previous value
	run_tt_gameplay(); // enable limited gameplay elements in tiled terrain mode
	draw_game_elements(timer1);
	draw_teleporters();
	if (camera_surf_collide && !no_tt_footsteps) {play_camera_footstep_sound();}
	if (change_near_far_clip) {check_zoom();} // reset perspective (may be unnecessary since will be reset on the next frame)
	check_xy_offsets();
	init_x = 0;
	if (TIMETEST) PRINT_TIME("3.9");
	//timing_profiler_stats();
}

