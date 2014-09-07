// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "mesh.h"
#include "main.h"
#include "timetest.h"
#include "physics_objects.h"
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


bool mesh_invalidated(1), no_asteroid_dust(0), fog_enabled(0);
int iticks(0), time0(0), scrolling(0), dx_scroll(0), dy_scroll(0), timer_a(0);
unsigned reflection_tid(0), enabled_lights(0); // 8 bit flags
float fticks(0.0), tfticks(0.0), tstep(0.0), camera_shake(0.0), cur_fog_end(1.0);
upos_point_type cur_origin(all_zeros);
colorRGBA cur_fog_color(GRAY);


extern bool nop_frame, combined_gu, have_sun, clear_landscape_vbo, show_lightning, spraypaint_mode, enable_depth_clamp, enable_multisample;
extern unsigned inf_terrain_fire_mode;
extern int auto_time_adv, camera_flight, reset_timing, run_forward, window_width, window_height, voxel_editing;
extern int advanced, b2down, dynamic_mesh_scroll, spectate, animate2, used_objs, disable_inf_terrain, curr_window, DISABLE_WATER;
extern float TIMESTEP, cloud_cover, univ_sun_rad, atmosphere, vegetation, zmin, zbottom, ztop, ocean_wave_height, brightness;
extern float def_atmosphere, def_vegetation;
extern double camera_zh;
extern point mesh_origin, surface_pos, univ_sun_pos, orig_cdir, sun_pos, moon_pos;
extern vector3d total_wind;
extern colorRGBA sun_color, bkg_color;
extern water_params_t water_params;
extern vector<camera_filter> cfilters;

void check_xy_offsets();
void post_window_redisplay();
void display_universe();
void display_inf_terrain(float uw_depth);
void update_temperature(bool verbose);
void update_sound_loops();
bool indir_lighting_updated();
point get_universe_display_camera_pos();
colorRGBA get_inf_terrain_mod_color();
void draw_gl_console();


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
	if (fabs(camera_shake) < 0.1) camera_shake = 0.0;
	assert(!dist_less_than(eye, center, TOLERANCE));
	assert(up_vector != zero_vector);
	fgLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up_vector.x, up_vector.y, up_vector.z);
}


void apply_camera_offsets(point const &camera) {

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
	if (camera_mode == 1) pos += surface_pos;
	if (camera_mode != 1 && combined_gu) pos += get_camera_pos(); // universe is always centered around the camera
	return pos;
}


point get_moon_pos() {

	point pos(moon_pos);
	if (camera_mode == 1) pos += surface_pos;
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
			if (rad < sun_radius) {
				color = sun_color;	
			}
			else {
				float blend_const(1.0 - (rad - sun_radius)/(outer_radius - sun_radius)); // 1.0 = full sun, 0.0 = no sun
				blend_color(color, sun_color, color, blend_const*blend_const, 0);
			}
		}
	}
	if (combined_gu) return color;

	if (light_factor < 0.6 && line_intersect_sphere(p1, v12, get_moon_pos(), moon_radius)) {
		color = texture_color(MOON_TEX);
	}
	return color;
}


void draw_stuff(int draw_uw, int timer1) {

	if (draw_uw) {
		draw_bubbles();
	}
	else {
		draw_splashes();
		draw_snow();
		draw_trees();
		render_voxel_data(0);
		check_gl_error(20);
		if (TIMETEST) PRINT_TIME("O");
		render_models(0);
		check_gl_error(21);
		if (TIMETEST) PRINT_TIME("P");
		check_gl_error(22);
		if (TIMETEST) PRINT_TIME("Q");
		draw_scenery(0, 1);
		draw_teleporters();
		check_gl_error(23);
		if (TIMETEST) PRINT_TIME("R");
		draw_coll_surfaces(0, 1);
		check_gl_error(24);
		if (TIMETEST) PRINT_TIME("S");
		draw_cracks_and_decals();
		if (TIMETEST) PRINT_TIME("S2");
		draw_transparent_object_groups();
		check_gl_error(25);
	}
}


void log_location(point const &pos) {

	static std::ofstream out;
	static bool inited;

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
		cout << "FPS: " << framerate << "  loc: (";
		camera.print();
		cout << ") @ frame " << frame_counter << endl;
		log_location(camera);
		show_framerate = 0;
	}
	if (display_framerate) {
		static int fr_counter(0);
		static float fr2(0.0);
		if (fr_counter%5 == 0) fr2 = framerate;
		draw_framerate(fr2);
		++fr_counter;
	}
}


float get_framerate(int &timer_b) {

	static float fr_average(0.0);
	timer_b = GET_TIME_MS();

	if (timer_b > timer_a) { // skip zero time frames
		float const framerate(1000.0/float(timer_b - timer_a));
		timer_a = timer_b;
		//return framerate;
		float const NUM_AVG = 5; // average over several frames
		fr_average = ((fr_average == 0.0) ? framerate : ((NUM_AVG - 1)*fr_average + framerate)/NUM_AVG);
	}
	return fr_average;
}


void final_draw(float framerate) {

	fog_enabled = 0;
	fgLoadIdentity();
	draw_camera_filters(cfilters);
	draw_frame_rate(framerate);
	show_other_messages();
}


void swap_buffers_and_redraw() {

	draw_gl_console();
	glutSwapBuffers();
	if (animate) {post_window_redisplay();} // before glutSwapBuffers()?
}


void calc_bkg_color() {

	float const lfn(1.0 - 5.0*(light_factor - 0.4));

	if (!have_sun || light_factor <= 0.4) {
		bkg_color = BACKGROUND_NIGHT;
	}
	else if (light_factor >= 0.6) {
		bkg_color = BACKGROUND_DAY;
	}
	else {
		blend_color(bkg_color, BACKGROUND_NIGHT, BACKGROUND_DAY, lfn, 1);
	}
	if (is_cloudy) {
		colorRGBA const orig_bkgc(bkg_color);
		blend_color(bkg_color, bkg_color, GRAY, 0.5, 1);
		UNROLL_3X(bkg_color[i_] = min(bkg_color[i_], orig_bkgc[i_]);) // can't make it brighter
	}
	if (atmosphere < 1.0) {blend_color(bkg_color, bkg_color, BACKGROUND_NIGHT, atmosphere, 0);}
}


float get_lf_scale(float lf) {
	return CLIP_TO_01(5.0f*(lf - 0.4f));
}

float get_moon_light_factor() {
	return fabs(moon_rot/PI - 1.0);
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


void config_bkg_color_and_clear(float depth, bool no_fog) {

	calc_bkg_color();
	glClearColor_rgba((!no_fog && show_fog) ? GRAY : bkg_color);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); // Clear the background
}


void reset_planet_defaults() {

	have_sun   = 1;
	atmosphere = def_atmosphere;
	vegetation = def_vegetation;
}


float get_light_pos_scale() {
	return ((world_mode == WMODE_INF_TERRAIN) ? 10.0 : 1.0); // hack: make the sun and moon far away in inf terrain mode 
}


void setup_lighting(float depth) {
	
	// background color code
	config_bkg_color_and_clear(depth, (world_mode == WMODE_INF_TERRAIN));

	// lighting code - RGB intensity for ambient and diffuse (specular is set elsewhere per object)
	float const mlf(get_moon_light_factor());
	colorRGBA ambient, diffuse;
	ambient[3] = diffuse[3] = 1.0;
	enabled_lights = 0;

	// Note: should this be set in universe lighting?
	colorRGBA ambient_c;
	blend_color(ambient_c, bkg_color, WHITE, 0.5, 1);

	for (unsigned i = 0; i < 3; ++i) {
		diffuse[i] = DSCALE*SUN_LT_C[i];
		ambient[i] = ASCALE*ambient_c[i];
	}
	for (unsigned i = 0; i < 3; ++i) {
		if (is_cloudy) {
			if (auto_time_adv) {
				diffuse[i] -= 0.06;
				ambient[i] -= 0.025;
			}
			else {
				diffuse[i] -= 0.15;
				ambient[i] -= 0.06;
			}
		}
		diffuse[i] -= 0.12*cloud_cover;
		ambient[i] -= 0.05*cloud_cover;
	}
	if (light_factor <= 0.4) { // moon
		calc_moon_atten(ambient, diffuse, mlf);
		sm_green_int = diffuse[1];
		set_colors_and_enable_light(1, ambient, diffuse);
	}
	else if (light_factor >= 0.6) { // sun
		for (unsigned i = 0; i < 3; ++i) { // add more brightness
			diffuse[i] += 0.2;
			ambient[i] += 0.1;
		}
		sm_green_int = diffuse[1];
		set_colors_and_enable_light(0, ambient, diffuse);
	}
	else { // sun and moon
		float const lfd(5.0*(light_factor - 0.4)), lfn(1.0 - lfd);

		for (unsigned i = 0; i < 3; ++i) { // should diffuse depend more on angle than ambient?
			diffuse[i] = (diffuse[i] + 0.2)*lfd;
			ambient[i] = (ambient[i] + 0.1)*lfd;
		}
		sm_green_int = lfd*diffuse[1];
		set_colors_and_enable_light(0, ambient, diffuse); // sun

		for (unsigned i = 0; i < 3; ++i) {
			diffuse[i] = (diffuse[i]/lfd - 0.2)*lfn;
			ambient[i] = (ambient[i]/lfd - 0.1)*lfn;
		}
		calc_moon_atten(ambient, diffuse, mlf);
		sm_green_int += lfn*diffuse[1];
		set_colors_and_enable_light(1, ambient, diffuse); // moon
	}

	// setup light position (after enabling lights)
	set_gl_light_pos(0, sun_pos *get_light_pos_scale(), LIGHT_W_VAL);
	set_gl_light_pos(1, moon_pos*get_light_pos_scale(), LIGHT_W_VAL);
	setup_gl_light_atten(0, 1.0, 0.0, 0.0); // reset attenuation to 1.0
}


void draw_sun_moon_stars() {

	if (light_factor <= 0.4) { // moon
		if (!is_cloudy) {gen_stars(1.0);}
		draw_moon();
	}
	else if (light_factor >= 0.6) { // sun
		draw_sun();
	}
	else { // sun and moon
		float const lfn(1.0 - 5.0*(light_factor - 0.4));
		if (!is_cloudy) {gen_stars(lfn);}
		draw_moon();
		draw_sun();
	}
}


void draw_universe_bkg(float depth, bool reflection_mode) {

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

	if (reflection_mode) { // setup mirror transform
		vector3d norm(0.0, 0.0, 1.0); // MZ
		rotate_vector3d_by_vr(new_sp, old_sp, norm);
		mirror_about_plane(norm, player_pos); // FIXME: doesn't use water_plane_z, so can't really be correct?
		//mirror_about_plane(norm, (player_pos + UNIV_NCLIP_SCALE*norm*(water_plane_z - camera.z)));
		set_player_dir(get_player_dir() - 2*norm*dot_product(get_player_dir(), norm)); // reflect the player view frustum dir
		set_univ_pdu();
	}

	// draw universe as background
	if (!reflection_mode) {config_bkg_color_and_clear(depth, 1);}
	point const camera_pos_orig(camera_pos);
	camera_pos = player_pos; // trick universe code into thinking the camera is at the player's ship
	stop_player_ship();
	if (TIMETEST) PRINT_TIME("0.1");
	bool const no_stars(is_cloudy || (atmosphere > 0.8 && light_factor >= 0.6));
	no_asteroid_dust = (reflection_mode || no_stars); // FIXME: should really pass this down (5 levels of function calls)
	draw_universe(1, 1, (no_stars ? 2 : 0)); // could clip by horizon?
	no_asteroid_dust = 0;
	if (TIMETEST) PRINT_TIME("0.2");
	camera_pos = camera_pos_orig;
	fgPopMatrix(); // undo universe transform

	// setup sun light source
	float const sun_intensity(max(0.25f, min(4.0f, 1000.0f*univ_sun_rad/sun_pos.mag())));
	setup_current_system(sun_intensity, reflection_mode);
	set_gl_light_pos(0, sun_pos*get_light_pos_scale(), LIGHT_W_VAL);
	disable_light(1); // no moonlight (for now)
	if (!have_sun || light_factor < 0.5) {set_light_ds_color(0, BLACK);} // sun below horizon: no diffuse or specular
	check_zoom(); // reset perspective
	//light_factor = (PI + 2.0*asinf(dot_product(sun_pos.get_norm(), plus_z)))/TWO_PI; // shouldn't change much from previous light_factor

	// setup background and init for standard mesh draw
	if (light_factor > 0.4) { // translucent blue for atmosphere
		colorRGBA color(bkg_color);
		if (reflection_mode) {color = cur_fog_color;} // make the sky reflection in the background blend with the fog in the foreground
		color.alpha *= 0.75*atmosphere*min(1.0, (light_factor - 0.4)/0.2);
		vector<camera_filter> cfs;
		cfs.push_back(camera_filter(color, 1, -1, 0));
		draw_camera_filters(cfs);
	}
	do_look_at();
	glClear(GL_DEPTH_BUFFER_BIT);
	if (TIMETEST) PRINT_TIME("0.3");
}


void draw_game_elements(int timer1) {

	if (TIMETEST) PRINT_TIME("U");
	draw_camera_weapon(1);
	draw_projectile_effects();
	if (TIMETEST) PRINT_TIME("V");
	draw_smoke_and_fires();
	if (TIMETEST) PRINT_TIME("W");
	draw_scheduled_weapons();
}


void setup_basic_fog() {

	if (!show_fog) return;
	setup_linear_fog(GRAY, 2.5*Z_SCENE_SIZE);
	fog_enabled = 1;
}


void add_uw_light_color_comp(int light, point const &lpos, float weight, colorRGBA &color) {

	// check for in shadow? what about tiled terrain?
	weight *= 0.5 + 0.5*max(0.0f, lpos.z/lpos.mag()); // vertical component (which penetrates water)
	UNROLL_3X(color[i_] += weight;)
}


void atten_uw_fog_color(colorRGBA &color, float depth) {

	colorRGBA light_color(BLACK);
	float const lf(get_lf_scale(light_factor));
	if (lf > 0.0) add_uw_light_color_comp(0, sun_pos,  lf,           light_color);
	if (lf < 1.0) add_uw_light_color_comp(1, moon_pos, 0.5*(1.0-lf), light_color);
	if (is_cloudy) light_color *= 0.5;
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
		atten_uw_fog_color(fog_color, 2.0*water_params.alpha*(water_plane_z - camera_z)); // more opaque = effectively deeper
		fog_dist = (0.3 + 1.5*Z_SCENE_SIZE*(camera_z - zmin2)/max(1.0E-3f, (water_plane_z - zmin2))) * max(0.1, (1.5 - water_params.alpha));
	}
	else {
		get_avg_sky_color(fog_color);
		colorRGBA const cloud_color(get_cloud_color(), 1.0); // alpha = 1.0
		blend_color(fog_color, cloud_color, bkg_color, 0.375, 1); // weighted more towards bkg_color
		fog_dist = get_inf_terrain_fog_dist();
	}
	setup_linear_fog(fog_color, fog_dist); // under water/ice
	fog_enabled = 1;
}


void scroll_scene() {

	RESET_TIME;
	point const camera(get_camera_pos());
	cout << "Shifting Scene..." << endl;
	camera_change = 1;
	scrolling     = 1;
	dx_scroll     = int(camera.x*DX_VAL_INV);
	dy_scroll     = int(camera.y*DY_VAL_INV);
	vector3d const vd(-DX_VAL*dx_scroll, -DY_VAL*dy_scroll, 0.0);
	surface_pos  += vd;
	xoff2        += dx_scroll;
	yoff2        += dy_scroll;
	shift_all_objs(vd);
	//PRINT_TIME("*** Top Level: Shift All Objects");
	reset_shadows(SHADOWED_ALL);
	gen_scene(1, 1, 1, 0, 0);
	//PRINT_TIME("*** Top Level: Gen Scene");
	regen_lightmap(); // not shiftable
	if (display_mode & 0x04) {water_plane_z = get_water_z_height();}
	update_temperature(0);
	recreated = 1;
	scrolling = 0;
	clear_landscape_vbo = 1;
	PRINT_TIME("*** Top Level: Final");
}


float get_ocean_wave_height() {

	if (!(display_mode & 0x0100)) return 0.0;
	static float time(0.0);
	if (animate2 && temperature > W_FREEZE_POINT) time += fticks;
	return ocean_wave_height*sin(1.0*time/TICKS_PER_SECOND); // add small waves
}


void draw_sun_flare() {

	//RESET_TIME;
	point const sun_pos(get_sun_pos());

	if (have_sun && light_factor >= 0.4 && sphere_in_camera_view(sun_pos, 4.0*sun_radius, 0)) { // use larger radius to include the flare/halo
		point const viewer(get_camera_pos());
		float intensity(1.0);

		if (world_mode == WMODE_GROUND) {
			unsigned const npts = 16;
			static point pts[npts];
			static bool pts_valid(0);
			unsigned nvis(0);
		
			for (unsigned i = 0; i < npts; ++i) {
				int index; // unused
				if (!pts_valid) {pts[i] = signed_rand_vector_norm();}
				point const pos(sun_pos + pts[i]*sun_radius);
				if (coll_pt_vis_test(pos, viewer, 0.0, index, camera_coll_id, 0, 1) && (!(display_mode & 0x01) || !line_intersect_mesh(pos, viewer, 0))) {++nvis;}
			}
			pts_valid = 1;
			if (nvis == 0) return;
			intensity = 0.1 + 0.9*float(nvis)/float(npts);
		}
		else if (world_mode == WMODE_INF_TERRAIN) {
			if (sun_pos.z < zmin) return; // sun below the mesh
			if (viewer.z < water_plane_z) {intensity = CLIP_TO_01(1.0f - 1.0f*(water_plane_z - viewer.z));} // attenuate sun flare when underwater
		}
		DoFlares(viewer, camera_origin, sun_pos, 1.0, (combined_gu ? 15.0*univ_sun_rad : 1.0), intensity);
	}
	//PRINT_TIME("Query + Flare");
}


void set_multisample(bool enable) {

	if (!enable_multisample) return;
	if (enable) {glEnable(GL_MULTISAMPLE);} else {glDisable(GL_MULTISAMPLE);}
}


// The display function. It is called whenever the window needs
// redrawing (ie: overlapping window moves, resize, maximize)
// display() is also called every so many milliseconds to provide a decent framerate
void display(void) {

	check_gl_error(0);
	if (glutGetWindow() != curr_window) return; // only process the current window

	if (nop_frame) { // force display sync after enter/leave game mode (or something like that)
		nop_frame = 0;
		return;
	}
	RESET_TIME;
	static int init(0), frame_index(0), time_index(0), global_time(0), tticks(0);
	static point old_spos(0.0, 0.0, 0.0);
	proc_kbd_events();

	if (!init) {
		init   = 1;
		fticks = 1.0;
		time0  = timer1;
	}
	else if (animate && !DETERMINISTIC_TIME) {
		float ftick(0.0);
		static float carry(0.0);
		float const time_delta((TICKS_PER_SECOND*(timer1 - time0))/1000.0);

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
			carry   = ftick - (float)iticks;
		}
		fticks = max(TOLERANCE, 0.9f*fticks + 0.1f*(ftick - carry)); // slow averaging filter
		assert(fticks >  0.0 && fticks < 1.0E12);
		assert(iticks >= 0   && iticks < 1000000000);
	}
	else {
		fticks = 1.0;
		iticks = 1;
	}
	tstep        = TIMESTEP*fticks;
	reset_timing = 0;
	check_gl_error(1);
	set_fill_mode();

	if (map_mode && world_mode != WMODE_UNIVERSE) {
		draw_overhead_map();

		if (world_mode == WMODE_INF_TERRAIN) { // map mode infinite terrain
			camera_origin = surface_pos;
			apply_camera_offsets(get_camera_pos());
		}
		swap_buffers_and_redraw();
		check_xy_offsets();

		if (world_mode == WMODE_GROUND) {
			process_groups(); // ???
			if (game_mode) update_game_frame(); // ???
		}
		return;
	}
	displayed = 1;
	up_vector.assign(0.0, sinf(up_theta), camera_y*cosf(up_theta));
	setup_sphere_vbos();
	check_gl_error(2);
	update_sound_loops();
	set_multisample(1);
	if (enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);} else {glDisable(GL_DEPTH_CLAMP);}

	if (world_mode == WMODE_UNIVERSE) {
		display_universe(); // infinite universe
	}
	else {
		if (!pause_frame) uevent_advance_frame();
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
		int timer_b;
		float const framerate(get_framerate(timer_b));
		if (global_time == 0) global_time = timer_b;
		
		if (show_framerate == 2) {
			cout << used_objs << " objects, time = " << (timer_b - global_time) << endl;
			cout << "Elapsed frames = " << (frame_counter - frame_index) << ", elapsed time = " << (timer_b - time_index)
				 << ", avg framerate = " << 1000.0*float(frame_counter - frame_index)/float(timer_b - time_index) << endl;
			frame_index    = frame_counter;
			time_index     = timer_b;
			show_framerate = 0;
		}
		if (show_framerate == 1) {
			timer_a = GET_TIME_MS();
			show_framerate = 2;
		}
		if (world_mode == WMODE_GROUND) {process_platforms();} // must be before camera code
		if (world_mode == WMODE_INF_TERRAIN) {camera_mode = 1;} // force to ground/walking mode

		// camera position code
		auto_advance_camera();

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
			if (c_radius >= C_RADIUS0) temp_c_radius = c_radius;
			c_radius      = CR_SCALE*C_RADIUS0;
			camera_origin = cpos2;
			if (camera_mode == 1 && !spectate) camera_origin.z += camera_zh;
		}
		else {
			if (temp_c_radius >= C_RADIUS0) {
				c_radius      = temp_c_radius;
				temp_c_radius = 0.0;
			}
		}
		update_cpos();
		point const camera(get_camera_pos());
		float depth;
		underwater = check_underwater(CAMERA_ID, depth);
		auto_advance_time();
		if (animate2) {total_wind += wind*fticks;}
		check_gl_error(3);
		if (TIMETEST) PRINT_TIME("\n\nA");
		
		if (!combined_gu) {
			do_look_at();
			sun_color = SUN_LT_C;
			apply_red_sky(sun_color);
			reset_planet_defaults();
			setup_lighting(depth);
			check_gl_error(4);
			if (TIMETEST) PRINT_TIME("B");
		}
		if (world_mode == WMODE_INF_TERRAIN) { // infinite terrain mode
			display_inf_terrain(depth);
		}
		else { // finite terrain mode
			if (combined_gu) { // light from current system's star
				draw_universe_bkg(depth, 0); // infinite universe as background
			}
			if (TIMETEST) PRINT_TIME("C");

			if (mesh_invalidated) {
				gen_mesh_bsp_tree();
				mesh_invalidated = 0;
			}
			// draw background
			if (!combined_gu) {draw_sun_moon_stars();}
			if (show_fog || underwater) {fog_enabled = 1;}
			if (!show_lightning) {l_strike.enabled = 0;}
			compute_brightness();
			if (!combined_gu) {draw_earth();}
			draw_sky(0);
			draw_puffy_clouds(0);
			check_gl_error(5);
			if (TIMETEST) PRINT_TIME("D");

			// run physics and collision detection
			reset_shadows(DYNAMIC_SHADOW);
			process_groups();
			check_gl_error(12);
			if (TIMETEST) PRINT_TIME("E");
			if (b2down) fire_weapon();
			update_weapon_cobjs(); // and update cblade
			check_gl_error(6);
			if (TIMETEST) PRINT_TIME("F");

			// create shadow map
			if (!camera_view) {camera_shadow(camera);}
			create_shadow_map(); // where should this go?
			if (TIMETEST) PRINT_TIME("G");

			proc_voxel_updates();

			// send data to GPU
			setup_object_render_data();
			check_gl_error(101);

			if (underwater) {
				colorRGBA fog_color(((temperature <= W_FREEZE_POINT) ? ICE_C : WATER_C), 1.0); // under ice/water, alpha = 1.0
				select_liquid_color(fog_color, camera);
				atten_uw_fog_color(fog_color, depth);
				float const fog_dist(0.2 + (0.25 + 0.75*fog_color.B)*(1.5*Z_SCENE_SIZE)*(camera.z - zmin)/((camera.z + depth) - zmin));
				setup_linear_fog(fog_color, fog_dist);
			}

			// draw the scene
			draw_camera_weapon(0);
			if (TIMETEST) PRINT_TIME("H");

			draw_coll_surfaces(1, 0);
			if (TIMETEST) PRINT_TIME("I");
			
			if (display_mode & 0x01) {display_mesh();} // draw mesh
			check_gl_error(7);
			if (TIMETEST) PRINT_TIME("J");

			draw_grass();
			draw_scenery(1, 0);
			if (TIMETEST) PRINT_TIME("K");

			draw_solid_object_groups();
			check_gl_error(8);
			if (TIMETEST) PRINT_TIME("L");

			draw_stuff(!underwater, timer1);

			if (show_lightning) { // after the water?
				l_strike.gen();

				if (l_strike.enabled == 1 && animate2) {
					if ((rand()&1) == 0) gen_smoke(l_strike.end);
					if ((rand()&7) == 0) gen_fire(l_strike.end, 1.0, NO_SOURCE);
				}
			}
			if (TIMETEST) PRINT_TIME("M");
			if (display_mode & 0x04) {draw_water();} // must be after process_groups()
			check_gl_error(9);
			if (TIMETEST) PRINT_TIME("N");
			draw_stuff(underwater, timer1);
			if (TIMETEST) PRINT_TIME("T");
			draw_game_elements(timer1);
			setup_basic_fog();
			draw_sky(1);
			draw_puffy_clouds(1);
			draw_sun_flare();
			check_gl_error(10);
			//draw_scene_bounds_and_light_frustum(get_light_pos()); // TESTING
		} // WMODE_GROUND
		check_gl_error(13);
		check_zoom();
		final_draw(framerate);
		purge_coll_freed(0); // optional
		camera_flight = 0;

		if (game_mode) {
			update_game_frame();
			show_user_stats();
			show_blood_on_camera();
			show_crosshair(WHITE, do_zoom);
		}
		else if (world_mode == WMODE_INF_TERRAIN && inf_terrain_fire_mode) {
			show_crosshair(get_inf_terrain_mod_color(), do_zoom);
		}
		else if (world_mode == WMODE_GROUND && voxel_editing) {
			show_crosshair(((voxel_editing == 2) ? RED : GREEN), do_zoom);
		}
		else if (spraypaint_mode) {
			draw_spraypaint_crosshair();
		}
		if (world_mode == WMODE_INF_TERRAIN || (world_mode == WMODE_GROUND && !game_mode && camera_mode == 1)) {
			draw_compass_and_alt();
		}
		if (indir_lighting_updated()) {
			draw_text(PURPLE, 0.007*(float)window_width/(float)window_height, -0.009, -0.02, "Lighting Updating");
		}
		if (TIMETEST) PRINT_TIME("X");

		if (dynamic_mesh_scroll && world_mode == WMODE_GROUND && camera_mode == 1 && !camera_view) {
			float const cdist(max(fabs(camera.x/X_SCENE_SIZE), fabs(camera.y/Y_SCENE_SIZE)));
			if (cdist > REL_SCROLL_DIST) scroll_scene();
		}
	} // not universe mode
	draw_enabled_ui_menus();
	swap_buffers_and_redraw();
	check_gl_error(11);
	if (TIMETEST) PRINT_TIME("Y");
}


void display_universe() { // infinite universe

	int timer_b;
	float framerate;
	static int init(0);
	RESET_TIME;

	if (!init || init_x) {
		init     = 1;
		init_x   = 0;
		show_fog = 0;
		update_cpos();
	}
	camera_view = 0;
	framerate   = get_framerate(timer_b);
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
	draw_blasts();
	if (TIMETEST) PRINT_TIME("Draw Blasts");
	final_draw(framerate);
	show_crosshair(WHITE, do_zoom);
	draw_universe_stats();
	camera_surf_collide = last_csc;
	check_gl_error(33);
}


void draw_transparent(bool above_water) {

	if (above_water) {
		draw_transparent_object_groups();
	}
	else {
		draw_bubbles();
	}
}


void apply_z_mirror(float zval) {

	fgTranslate(0.0, 0.0, 2*zval); // translate to zval and back
	fgScale(1.0, 1.0, -1.0); // scale in z
	//mirror_about_plane(plus_z, point(0.0, 0.0, zval));
}


// render scene reflection to texture
void create_reflection_texture(unsigned tid, unsigned xsize, unsigned ysize, float terrain_zmin) {

	//RESET_TIME;
	// setup reflected camera frustum
	pos_dir_up const old_camera_pdu(camera_pdu); // reflect camera frustum used for VFC
	camera_pdu.pos.z = 2*water_plane_z - camera_pdu.pos.z;
	camera_pdu.dir.z = -camera_pdu.dir.z; // mirror
	camera_pdu.upv_  = -camera_pdu.upv_;
	camera_pdu.orthogonalize_up_dir();
	pos_dir_up const refl_camera_pdu(camera_pdu);

	pre_draw_tiled_terrain();

	// setup viewport and projection matrix
	glViewport(0, 0, xsize, ysize);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	if (combined_gu && !is_cloudy) {draw_universe_bkg(0.0, 1);} // infinite universe as background
	fgMatrixMode(FG_PROJECTION);
	fgPushMatrix();
	set_perspective(PERSP_ANGLE, 1.0);
	do_look_at();

	// setup mirror transform
	fgMatrixMode(FG_MODELVIEW);
	fgPushMatrix();
	apply_z_mirror(water_plane_z);
	camera_pdu = refl_camera_pdu; // reset reflected PDU

	// draw partial scene
	if (!combined_gu) {
		draw_sun_moon_stars();
		draw_sun_flare();
	}
	draw_cloud_planes(terrain_zmin, 1, 1, 0); // slower but a nice effect
	if (show_lightning) {draw_tiled_terrain_lightning(1);}
	if (get_camera_pos().z <= get_tt_cloud_level()) {draw_tiled_terrain(1);} // camera is below the clouds
	fgPopMatrix(); // end mirror transform

	// render reflection to texture
	bind_2d_texture(tid);
	glReadBuffer(GL_BACK);
	// glCopyTexSubImage2D copies the frame buffer to the bound texture
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, xsize, ysize);

	// reset state
	camera_pdu = old_camera_pdu; // restore camera_pdu
	fgMatrixMode(FG_PROJECTION);
	fgPopMatrix();
	fgMatrixMode(FG_MODELVIEW);
	set_standard_viewport();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//PRINT_TIME("Create Reflection Texture");
}


unsigned create_reflection(float terrain_zmin) {

	if (display_mode & 0x20) return 0; // reflections not enabled
	static unsigned last_xsize(0), last_ysize(0);
	unsigned const xsize(window_width/2), ysize(window_height/2);

	if (last_xsize != xsize || last_ysize != ysize) {
		free_texture(reflection_tid);
		last_xsize = xsize;
		last_ysize = ysize;
	}
	if (!reflection_tid) {
		setup_texture(reflection_tid, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, xsize, ysize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	}
	assert(glIsTexture(reflection_tid));
	create_reflection_texture(reflection_tid, xsize, ysize, terrain_zmin);
	check_gl_error(999);
	return reflection_tid;
}


void display_inf_terrain(float uw_depth) { // infinite terrain mode (Note: uses light params from ground mode)

	static int init_xx(1);
	RESET_TIME;

	if (init_x || init_xx) {
		init_xx  = 0;
		show_fog = 1;
		c_radius = 0.001;
		update_cpos();
		surface_pos         = all_zeros;
		camera_surf_collide = 1;
	}
	camera_view = 0;
	if (camera_surf_collide) {check_player_tiled_terrain_collision();}
	update_temperature(0);
	point const camera(get_camera_pos());
	apply_camera_offsets(camera);
	compute_brightness();
	set_global_state();
	if (b2down) fire_weapon();

	bool const water_enabled((display_mode & 0x04) && !DISABLE_WATER);
	water_plane_z = (water_enabled ? (get_water_z_height() + get_ocean_wave_height()) : -10*FAR_CLIP);
	camera_mode   = 1; // walking on ground
	float min_camera_dist(0.0);
	float const terrain_zmin(update_tiled_terrain(min_camera_dist));
	bool const change_near_far_clip(!camera_surf_collide && min_camera_dist > 0.0);
	bool const draw_water(water_enabled && water_plane_z >= terrain_zmin);
	if (show_fog || underwater) {set_inf_terrain_fog(underwater, terrain_zmin);}
	unsigned reflection_tid(0);
	if (draw_water && !underwater) {reflection_tid = create_reflection(terrain_zmin);}

	if (combined_gu) {
		draw_universe_bkg(uw_depth, 0); // infinite universe as background
		check_gl_error(4);
	}
	else {
		config_bkg_color_and_clear(uw_depth, 1);
		draw_sun_moon_stars();
	}
	if (change_near_far_clip) {
		float const near_clip(NEAR_CLIP + 0.01*min_camera_dist);
		float const far_clip(get_tt_fog_based_far_clip(min_camera_dist));
		set_perspective_near_far(near_clip, far_clip);
		do_look_at(); // clear depth buffer?
		camera_pdu.near_ = near_clip; // override camera frustum near/far clip so that VFC will be correct
		camera_pdu.far_  = far_clip;
	}
	bool const camera_above_clouds(camera.z > get_tt_cloud_level());
	draw_cloud_planes(terrain_zmin, 0, !camera_above_clouds, 1); // these two lines could go in either order
	draw_sun_flare();
	if (TIMETEST) PRINT_TIME("3.2");
	if (show_lightning) {draw_tiled_terrain_lightning(0);}
	pre_draw_tiled_terrain();
	if (TIMETEST) PRINT_TIME("3.26");
	draw_tiled_terrain(0);
	if (TIMETEST) PRINT_TIME("3.3");
	//if (underwater ) {draw_tiled_terrain_precipitation();}
	if (draw_water ) {draw_water_plane(water_plane_z, reflection_tid);}
	if (!underwater) {draw_tiled_terrain_precipitation();}
	draw_cloud_planes(terrain_zmin, 0, camera_above_clouds, 0);
	if (change_near_far_clip) {check_zoom();} // reset perspective (may be unnecessary since will be reset on the next frame)
	check_xy_offsets();
	init_x = 0;
	if (TIMETEST) PRINT_TIME("3.9");
}

