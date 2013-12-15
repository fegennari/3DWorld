// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#include "mesh.h"
#include "main.h"
#include "timetest.h"
#include "physics_objects.h"
#include "model3d.h"
#include <fstream>


/* GL_LIGHT<N>:
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


bool mesh_invalidated(1), no_asteroid_dust(0);
int iticks(0), time0(0), scrolling(0), dx_scroll(0), dy_scroll(0), timer_a(0);
unsigned reflection_tid(0);
float fticks(0.0), tfticks(0.0), tstep(0.0), camera_shake(0.0);
upos_point_type cur_origin(all_zeros);


extern bool nop_frame, combined_gu, have_sun, clear_landscape_vbo, show_lightning, spraypaint_mode;
extern unsigned inf_terrain_fire_mode;
extern int auto_time_adv, camera_flight, reset_timing, enable_fsource, run_forward, window_width, window_height, voxel_editing;
extern int advanced, b2down, dynamic_mesh_scroll, spectate, animate2, used_objs, disable_inf_terrain, curr_window, DISABLE_WATER;
extern float TIMESTEP, cloud_cover, univ_sun_rad, atmosphere, vegetation, zmin, zbottom, ztop, ocean_wave_height;
extern double camera_zh;
extern point mesh_origin, flow_source, surface_pos, univ_sun_pos, orig_cdir, sun_pos, moon_pos;
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
void calc_cur_ambient_diffuse();
bool indir_lighting_updated();
colorRGBA get_tt_water_color();


void glClearColor_rgba(const colorRGBA &color) {

	glClearColor(color.R, color.G, color.B, color.A);
}


void set_fsource_light() {

	if (!enable_fsource) return;
	float const diffuse[4] = {1.0,  1.0,  1.0,  1.0};
	float const ambient[4] = {0.25, 0.25, 0.25, 1.0};
	int const light(GL_LIGHT2);
	set_colors_and_enable_light(light, ambient, diffuse);
	set_gl_light_pos(light, flow_source, 1.0);
	glLightf(light, GL_CONSTANT_ATTENUATION,  0.6);
	glLightf(light, GL_LINEAR_ATTENUATION,    1.2);
	glLightf(light, GL_QUADRATIC_ATTENUATION, 0.6);
}


void calc_moon_atten(float *ambient, float *diffuse, float mlf) {

	if (mlf < 0.6) {
		float const moon_atten((mlf < 0.5) ? 0.0 : 10.0*(mlf - 0.5));
		for (unsigned i = 0; i < 3; ++i) diffuse[i] *= moon_atten;
	}
	for (unsigned i = 0; i < 3; ++i) {
		ambient[i] *= 0.5;
		diffuse[i] *= 0.5;
	}
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
	gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up_vector.x, up_vector.y, up_vector.z);
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
		draw_snow();
		draw_trees();
		check_gl_error(20);
		if (TIMETEST) PRINT_TIME("O");
		draw_hmv();
		draw_scenery(0, 1);
		draw_teleporters();
		check_gl_error(21);
		if (TIMETEST) PRINT_TIME("P");
		draw_coll_surfaces(0, 1);
		check_gl_error(22);
		if (TIMETEST) PRINT_TIME("Q");
		render_models(0);
		check_gl_error(24);
		if (TIMETEST) PRINT_TIME("R");
		render_voxel_data(0);
		check_gl_error(25);
		if (TIMETEST) PRINT_TIME("S");
		draw_cracks_and_decals();
		if (TIMETEST) PRINT_TIME("S2");
		draw_transparent_object_groups();
		check_gl_error(26);
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
		point const camera(get_camera_pos());
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

	glDisable(GL_FOG);
	glLoadIdentity();
	draw_camera_filters(cfilters);
	draw_frame_rate(framerate);
	show_other_messages();
}


void swap_buffers_and_redraw() {

	glutSwapBuffers();
	if (animate) post_window_redisplay();
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

	for (unsigned i = 0; i < 4; ++i) {
		color[i] *= sun_color[i];
	}
	float cmult(have_sun ? get_lf_scale(light_factor) : 0.0); // light from sun

	if (!combined_gu && light_factor < 0.6) { // light from moon
		float const lfs(get_lf_scale(get_moon_light_factor()));
		cmult += 0.5*CLIP_TO_01(5.0f*(0.6f - light_factor))*lfs;
	}
	color *= (0.8*cmult + 0.2);
}


colorRGBA set_lighted_fog_color(colorRGBA color) {

	add_sun_effect(color);
	glFogfv(GL_FOG_COLOR, (float *)&color);
	return color;
}


void auto_advance_camera() {

	if (run_forward && !advanced) advance_camera(MOVE_FRONT);
	advanced = 0;
}


void config_bkg_color_and_clear(bool underwater, float depth, bool no_fog) {

	calc_bkg_color();
	glClearColor_rgba((!no_fog && show_fog) ? GRAY : bkg_color);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT /*| GL_ACCUM_BUFFER_BIT*/); // Clear the background
}


void reset_planet_defaults() {

	have_sun   = 1;
	atmosphere = 1.0;
	vegetation = 1.0;
}


void set_uniform_atten_lighting(int light) {

	glLightf(light, GL_CONSTANT_ATTENUATION,  1.0);
	glLightf(light, GL_LINEAR_ATTENUATION,    0.0);
	glLightf(light, GL_QUADRATIC_ATTENUATION, 0.0);
}


void setup_lighting(bool underwater, float depth) {
	
	// background color code
	config_bkg_color_and_clear(underwater, depth, (world_mode == WMODE_INF_TERRAIN));
	
	// setup light position
	float const light_pos_scale((world_mode == WMODE_INF_TERRAIN) ? 10.0 : 1.0); // hack: make the sun and moon far away in inf terrain mode 
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	set_gl_light_pos(GL_LIGHT0, sun_pos*light_pos_scale,  LIGHT_W_VAL);
	set_gl_light_pos(GL_LIGHT1, moon_pos*light_pos_scale, LIGHT_W_VAL);
	set_uniform_atten_lighting(GL_LIGHT0); // reset attenuation to 1.0

	// lighting code - RGB intensity for ambient and diffuse (specular is set elsewhere per object)
	set_fsource_light();
	float const mlf(get_moon_light_factor());
	float ambient[4], diffuse[4];
	ambient[3] = diffuse[3] = 1.0;

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
		set_colors_and_enable_light(GL_LIGHT1, ambient, diffuse);
	}
	else if (light_factor >= 0.6) { // sun
		for (unsigned i = 0; i < 3; ++i) { // add more brightness
			diffuse[i] += 0.2;
			ambient[i] += 0.1;
		}
		sm_green_int = diffuse[1];
		set_colors_and_enable_light(GL_LIGHT0, ambient, diffuse);
	}
	else { // sun and moon
		float const lfd(5.0*(light_factor - 0.4)), lfn(1.0 - lfd);

		for (unsigned i = 0; i < 3; ++i) { // should diffuse depend more on angle than ambient?
			diffuse[i] = (diffuse[i] + 0.2)*lfd;
			ambient[i] = (ambient[i] + 0.1)*lfd;
		}
		sm_green_int = lfd*diffuse[1];
		set_colors_and_enable_light(GL_LIGHT0, ambient, diffuse); // sun

		for (unsigned i = 0; i < 3; ++i) {
			diffuse[i] = (diffuse[i]/lfd - 0.2)*lfn;
			ambient[i] = (ambient[i]/lfd - 0.1)*lfn;
		}
		calc_moon_atten(ambient, diffuse, mlf);
		sm_green_int += lfn*diffuse[1];
		set_colors_and_enable_light(GL_LIGHT1, ambient, diffuse); // moon
	}
}


void draw_sun_moon_stars() {

	if (light_factor <= 0.4) { // moon
		if (!is_cloudy) {gen_stars(1.0, island);}
		draw_moon();
	}
	else if (light_factor >= 0.6) { // sun
		draw_sun();
	}
	else { // sun and moon
		float const lfn(1.0 - 5.0*(light_factor - 0.4));
		if (!is_cloudy) {gen_stars(lfn, island);}
		draw_moon();
		draw_sun();
	}
}


void draw_universe_bkg(bool underwater, float depth, bool reflection_mode) {

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
	glPushMatrix();
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
	if (!reflection_mode) {config_bkg_color_and_clear(underwater, depth, 1);}
	point const camera_pos_orig(camera_pos);
	camera_pos = player_pos; // trick universe code into thinking the camera is at the player's ship
	stop_player_ship();
	if (TIMETEST) PRINT_TIME("0.1");
	bool const no_stars(is_cloudy || (atmosphere > 0.8 && light_factor >= 0.6));
	int const fog_enabled(glIsEnabled(GL_FOG));
	if (fog_enabled) {glDisable(GL_FOG);}
	no_asteroid_dust = (reflection_mode || no_stars); // FIXME: should really pass this down (5 levels of function calls)
	draw_universe(1, 1, (no_stars ? 2 : 0)); // could clip by horizon?
	no_asteroid_dust = 0;
	if (fog_enabled) {glEnable(GL_FOG);}
	if (TIMETEST) PRINT_TIME("0.2");
	camera_pos = camera_pos_orig;
	glPopMatrix(); // undo universe transform

	// setup sun light source
	setup_current_system();
	glPushMatrix();
	translate_to(get_camera_pos());
	set_gl_light_pos(GL_LIGHT0, sun_pos, LIGHT_W_VAL);
	glPopMatrix();
	set_light_atten(GL_LIGHT0, max(0.25f, min(4.0f, 0.0007f*sun_pos.mag()/univ_sun_rad)));
	glDisable(GL_LIGHT1); // no moonlight (for now)
	
	if (!have_sun || light_factor < 0.5) { // sun below horizon
		glLightfv(GL_LIGHT0, GL_DIFFUSE, &BLACK.R); // no diffuse
	}
	check_zoom(); // reset perspective
	//light_factor = (PI + 2.0*asinf(dot_product(sun_pos.get_norm(), plus_z)))/TWO_PI; // shouldn't change much from previous light_factor

	// setup background and init for standard mesh draw
	if (light_factor > 0.4) { // translucent blue for atmosphere
		colorRGBA color(bkg_color);
		if (reflection_mode) {glGetFloatv(GL_FOG_COLOR, (float *)&color);} // make the sky reflection in the background blend with the fog in the foreground
		color.alpha *= 0.75*atmosphere*min(1.0, (light_factor - 0.4)/0.2);
		vector<camera_filter> cfs;
		cfs.push_back(camera_filter(color, 1, -1));
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
	set_lighted_fog_color(GRAY);
	glFogf(GL_FOG_END, 2.5*Z_SCENE_SIZE);
	glEnable(GL_FOG);
}


void add_uw_light_color_comp(int light, point const &lpos, float weight, colorRGBA &color) {

	// check for in shadow? what about tiled terrain?
	weight *= 0.5 + 0.5*max(0.0f, lpos.z/lpos.mag()); // vertical component (which penetrates water)
	UNROLL_3X(color[i_] += weight;)
}


void atten_uw_fog_color(colorRGBA &color, float depth) {

	colorRGBA light_color(BLACK);
	float const lf(get_lf_scale(light_factor));
	if (lf > 0.0) add_uw_light_color_comp(GL_LIGHT0, sun_pos,  lf,           light_color);
	if (lf < 1.0) add_uw_light_color_comp(GL_LIGHT1, moon_pos, 0.5*(1.0-lf), light_color);
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
		fog_color = get_tt_water_color();
		fog_color.alpha = 1.0;
		atten_uw_fog_color(fog_color, 2.0*water_params.alpha*(water_plane_z - camera_z)); // more opaque = effectively deeper
		fog_dist = max(0.1, (0.3 + 1.5*Z_SCENE_SIZE*(camera_z - zmin2)/max(1.0E-3f, (water_plane_z - zmin2))) * (1.5 - water_params.alpha)); // not sure why the max of 0.1 is required
	}
	else {
		get_avg_sky_color(fog_color);
		fog_dist = get_inf_terrain_fog_dist();
	}
	set_lighted_fog_color(fog_color); // under water/ice
	glFogf(GL_FOG_END, fog_dist);
	glFogf(GL_FOG_DENSITY, 0.05); // density isn't used, but it doesn't hurt to set it to around this value in case it is
	glEnable(GL_FOG);
}


void scroll_scene() {

	RESET_TIME;
	point const camera(get_camera_pos());
	cout << "Shifting Scene..." << endl;
	mesh_type     = 0; // don't scroll to an island
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
	//PRINT_TIME("*** Top Level: Shift All Objects + Reset Shadows");
	gen_scene(1, 1, 1, 0, 0);
	//PRINT_TIME("*** Top Level: Gen Scene");
	shift_lightmap(vd);
	if (display_mode & 0x04) water_plane_z = get_water_z_height();
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

	if (have_sun && light_factor >= 0.4 && sphere_in_camera_view(sun_pos, 4.0*sun_radius, 2)) { // use larger radius to include the flare/halo
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
				if (coll_pt_vis_test(pos, viewer, 0.0, index, camera_coll_id, 0, 1) && !line_intersect_mesh(pos, viewer, 0)) {++nvis;}
			}
			pts_valid = 1;
			if (nvis == 0) return;
			intensity = 0.1 + 0.9*float(nvis)/float(npts);
		}
		int const fog_enbaled(glIsEnabled(GL_FOG));
		glDisable(GL_FOG);
		DoFlares(viewer, camera_origin, sun_pos, 1.0, (combined_gu ? 15.0*univ_sun_rad : 1.0), intensity);
		if (fog_enbaled) glEnable(GL_FOG);
	}
	//PRINT_TIME("Query + Flare");
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
	set_lighted_sides(1);
	displayed = 1;
	up_vector.assign(0.0, sinf(up_theta), camera_y*cosf(up_theta));
	setup_sphere_vbos();
	check_gl_error(2);
	update_sound_loops();

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
			cout << used_objs << " objects + " << 2*(MESH_X_SIZE/resolution)*(MESH_Y_SIZE/resolution)
				 << " triangles, time = " << (timer_b - global_time) << endl;
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
			setup_lighting(underwater, depth);
			set_light_atten(GL_LIGHT0, 1.0); // reset light0
			check_gl_error(4);
			if (TIMETEST) PRINT_TIME("B");
		}
		if (world_mode == WMODE_INF_TERRAIN) { // infinite terrain mode
			display_inf_terrain(depth);
		}
		else { // finite terrain mode
			if (combined_gu) { // light from current system's star
				draw_universe_bkg(underwater, depth, 0); // infinite universe as background
			}
			if (TIMETEST) PRINT_TIME("C");

			if (mesh_invalidated) {
				gen_mesh_bsp_tree();
				mesh_invalidated = 0;
			}
			// draw background
			if (!combined_gu) {draw_sun_moon_stars();}
			if (show_fog || underwater) {glEnable(GL_FOG);}
			if (!show_lightning) {l_strike.enabled = 0;}
			compute_brightness();
			if (!combined_gu) {draw_earth();}
			draw_sky(0);
			draw_puffy_clouds(0);
			draw_env_other();
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
				bool const is_ice(temperature <= W_FREEZE_POINT && (!island || camera.z > ocean.z));
				colorRGBA fog_color(is_ice ? ICE_C : WATER_C); // under ice/water
				fog_color.alpha = 1.0;
				select_liquid_color(fog_color, camera);
				atten_uw_fog_color(fog_color, depth);
				set_lighted_fog_color(fog_color);
				float const fog_dist(0.2 + (0.25 + 0.75*fog_color.B)*(1.5*Z_SCENE_SIZE)*(camera.z - zmin)/((camera.z + depth) - zmin));
				glFogf(GL_FOG_END, fog_dist);
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
					if ((rand()&7) == 0) gen_fire(l_strike.end, 1.0, -2);
				}
			}
			if (TIMETEST) PRINT_TIME("M");
			if (display_mode & 0x04) {draw_water();} // must be after process_groups()
			if (display_mode & 0x01) {draw_ocean();}
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
			show_crosshair(((inf_terrain_fire_mode == 2) ? RED : ((inf_terrain_fire_mode == 3) ? BLUE : GREEN)), do_zoom);
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
			PURPLE.do_glColor();
			draw_text(0.007*(float)window_width/(float)window_height, -0.009, -0.02, "Lighting Updating");
		}
		if (TIMETEST) PRINT_TIME("X");

		if (dynamic_mesh_scroll && world_mode == WMODE_GROUND && camera_mode == 1 && !island && !camera_view) {
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
		mesh_type = INIT_MESH_TYPE;
		init      = 1;
		init_x    = 0;
		show_fog  = 0;
		update_cpos();
	}
	camera_view = 0;
	ocean_set   = 0;
	framerate   = get_framerate(timer_b);
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
	if (!pause_frame) uevent_advance_frame();
	do_look_at();
	if (b2down) fire_weapon();
	update_blasts();
	if (TIMETEST) PRINT_TIME("Process BRs");
	check_gl_error(31);
	draw_universe();
	check_gl_error(32);
	if (TIMETEST) PRINT_TIME("Draw Universe");
	draw_blasts();
	if (TIMETEST) PRINT_TIME("Draw Blasts");
	set_light_atten(GL_LIGHT0, 1.0); // reset light0
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

	glTranslatef(0.0, 0.0, 2*zval); // translate to zval and back
	glScalef(1.0, 1.0, -1.0); // scale in z
	//mirror_about_plane(plus_z, point(0.0, 0.0, zval));
}


// render scene reflection to texture
void create_reflection_texture(unsigned tid, unsigned xsize, unsigned ysize) {

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
	if (combined_gu && !is_cloudy) {draw_universe_bkg(underwater, 0.0, 1);} // infinite universe as background
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	set_perspective(PERSP_ANGLE, 1.0);
	do_look_at();

	// setup mirror transform
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	apply_z_mirror(water_plane_z);
	camera_pdu = refl_camera_pdu; // reset reflected PDU

	// draw partial scene
	if (!combined_gu) {
		draw_sun_moon_stars();
		draw_sun_flare();
	}
	draw_cloud_plane(zmin, 1); // slower but a nice effect
	if (show_lightning) {draw_tiled_terrain_lightning(1);}
	// setup above-water clip plane for mesh
	double const plane[4] = {0.0, 0.0, 1.0, -water_plane_z}; // water at z=-water_z (mirrored)
	glEnable(WATER_CLIP_PLANE);
	glClipPlane(WATER_CLIP_PLANE, plane);
	draw_tiled_terrain(1);
	glDisable(WATER_CLIP_PLANE);
	// could render more of the scene here
	glPopMatrix(); // end mirror transform

	// render reflection to texture
	bind_2d_texture(tid);
	glReadBuffer(GL_BACK);
	// glCopyTexSubImage2D copies the frame buffer to the bound texture
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, xsize, ysize);

	// reset state
	camera_pdu = old_camera_pdu; // restore camera_pdu
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	set_standard_viewport();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//PRINT_TIME("Create Reflection Texture");
}


unsigned create_reflection() {

	if (display_mode & 0x20) return 0; // reflections not enabled
	static unsigned last_xsize(0), last_ysize(0);
	unsigned const xsize(window_width/2), ysize(window_height/2);

	if (last_xsize != xsize || last_ysize != ysize) {
		free_texture(reflection_tid);
		last_xsize = xsize;
		last_ysize = ysize;
	}
	if (!reflection_tid) {
		setup_texture(reflection_tid, GL_MODULATE, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, xsize, ysize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	}
	assert(glIsTexture(reflection_tid));
	create_reflection_texture(reflection_tid, xsize, ysize);
	glDisable(GL_TEXTURE_2D);
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

	if (init_x && mesh_type != 0) {
		mesh_type = 0;
		gen_scene(1, 0, KEEP_MESH, 0, 0);
		recreated = 1;
	}
	bool const water_enabled((display_mode & 0x04) && !DISABLE_WATER);
	water_plane_z = (water_enabled ? (get_water_z_height() + get_ocean_wave_height()) : -10*FAR_CLIP);
	ocean.z       = water_plane_z;
	camera_mode   = 1; // walking on ground
	mesh_type     = 0;
	float const zmin2(update_tiled_terrain());
	bool const draw_water(water_enabled && water_plane_z >= zmin2);
	if (show_fog || underwater) {set_inf_terrain_fog(underwater, zmin2);}
	unsigned reflection_tid(0);
	if (draw_water && !underwater) reflection_tid = create_reflection();

	if (combined_gu) {
		draw_universe_bkg(underwater, uw_depth, 0); // infinite universe as background
		check_gl_error(4);
	}
	else {
		config_bkg_color_and_clear(underwater, uw_depth, 1);
		int const fog_enabled(glIsEnabled(GL_FOG));
		if (fog_enabled) {glDisable(GL_FOG);}
		draw_sun_moon_stars();
		if (fog_enabled) {glEnable(GL_FOG);}
	}
	draw_cloud_plane(zmin2, 0); // these two lines could go in either order
	draw_sun_flare();
	if (TIMETEST) PRINT_TIME("3.2");
	calc_cur_ambient_diffuse();
	if (TIMETEST) PRINT_TIME("3.25");
	if (show_lightning) {draw_tiled_terrain_lightning(0);}
	pre_draw_tiled_terrain();
	if (TIMETEST) PRINT_TIME("3.26");
	draw_tiled_terrain(0);
	if (TIMETEST) PRINT_TIME("3.3");
	if (underwater ) {draw_tiled_terrain_precipitation();}
	if (draw_water ) {draw_water_plane(water_plane_z, reflection_tid);}
	if (!underwater) {draw_tiled_terrain_precipitation();}
	check_xy_offsets();
	init_x = 0;
	if (TIMETEST) PRINT_TIME("3.9");
}

