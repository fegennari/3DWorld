// 3D World - OpenGL CS184 Computer Graphics Project
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
#include <set>

using namespace std;
typedef set<unsigned char>::iterator keyset_it;


int const DEF_SD           = 0; // shadow detail, max is 6
int const INIT_DMODE       = 0x010F;
int const P_MOTION_DEF     = 0;
int const KBD_HANDLER      = 1;
int const STARTING_INIT_X  = 0; // setting to 1 seems safer but less efficient
int const MIN_TIME_MS      = 100;

float const DEF_CRADIUS    = 10.0;
float const DEF_CTHETA     = -1.0;
float const DEF_CPHI       = 1.5;
float const DEF_UPTHETA    = 0.0;
float const DEF_CAMY       = 1.0;
float const MAP_SHIFT      = 8.0;
float const MAP_ZOOM       = 1.5;
float const D_TIMESTEP     = 1.5;
float const MOUSE_R_ADJ    = 0.01;   // mouse zoom radius adjustment
float const MOUSE_TRAN_ADJ = 0.01;   // mouse translate adjustment
float const MOUSE_ANG_ADJ  = PI/350; // mouse camera angle increment
float const MA_TOLERANCE   = 0.0001; // tolerance adjustment
float const CAMERA_AIR_CONT= 0.5;


char *defaults_file    = "defaults.txt";
char *dstate_file      = "state.txt";
char *dmesh_file       = "mesh.txt";
char *dcoll_obj_file   = "coll_objs/coll_objs.txt";
char *dmh_filename_raw = "heightmap.raw";
char *dmh_filename_bmp = "heightmap.bmp";
char *ddem_filename    = "dem/bot_fg.dem.txt";
char *ddem_raw_out     = "dem/bot_dem2.raw";
char *dship_def_file   = "ship_defs.txt";
char *state_file(dstate_file), *mesh_file(dmesh_file), *coll_obj_file(dcoll_obj_file);
char *mh_filename_raw(dmh_filename_raw), *mh_filename_bmp(dmh_filename_bmp), *dem_filename(ddem_filename);
char *dem_raw_out(ddem_raw_out), *ship_def_file(dship_def_file);
char *lighting_file(NULL), *lighting_file_l(NULL), *snow_file(NULL);


// Global Variables
bool nop_frame(0), combined_gu(0), underwater(0), kbd_text_mode(0), use_stencil_shadows(0), univ_stencil_shadows(1);
bool univ_planet_lod(0), draw_mesh_shader(1), show_lightning(0), disable_shaders(0), use_waypoints(0), group_back_face_cull(0);
int xoff(0), yoff(0), xoff2(0), yoff2(0), rand_gen_index(0), camera_change(1), camera_in_air(0), auto_time_adv(0);
int animate(1), animate2(1), begin_motion(0), draw_model(0), init_x(STARTING_INIT_X), fire_key(0), do_run(0);
int game_mode(0), map_mode(0), load_hmv(0), load_coll_objs(1), read_landscape(0), screen_reset(0), mesh_seed(0);
int display_framerate(1), init_resize(1), temp_change(0), mesh_type(INIT_MESH_TYPE), mt2(0), is_cloudy(0);
int star_init(0), recreated(1), stencil_shadow(0);
int displayed(0), min_time(0), resolution(1+(START_MODE==3)), res_old(1+(START_MODE!=3)), show_framerate(0);
int camera_view(0), camera_reset(1), camera_mode(0), camera_surf_collide(1), camera_coll_smooth(0);
int window_width(0), window_height(0), ww2(0), wh2(0), map_color(1); // window dimensions, etc.
int border_height(20), border_width(4), world_mode(START_MODE), display_mode(INIT_DMODE), do_read_mesh(0);
int last_mouse_x(0), last_mouse_y(0), m_button(0), mouse_state(1), maximized(0), verbose_mode(0), leaf_color_changed(0);
int shadow_detail(DEF_SD), do_zoom(0), disable_universe(0), disable_inf_terrain(0);
int num_trees(0), num_smileys(1), gmww = 640, gmwh = 480, srand_param(3), left_handed(0), mesh_scale_change(0);
int pause_frame(0), show_fog(0), spectate(0), b2down(0), free_for_all(0), teams(2), show_scores(0), universe_only(0);
int reset_timing(0), read_heightmap(0), default_ground_tex(-1), num_dodgeballs(1), INIT_DISABLE_WATER, ground_effects_level(2);
int enable_fsource(0), run_forward(0), advanced(0), passive_motion(P_MOTION_DEF), dynamic_mesh_scroll(0);
int read_light_file(0), write_light_file(0), read_light_file_l(0), write_light_file_l(0), read_snow_file(0), write_snow_file(0);
unsigned num_snowflakes(0), scenery_extent(0);
float water_plane_z(0.0), base_gravity(1.0), crater_size(1.0), tree_size(1.0), disabled_mesh_z(FAR_CLIP), vegetation(1.0);
float mesh_file_scale(1.0), mesh_file_tz(0.0), speed_mult(1.0), mesh_z_cutoff(-FAR_CLIP), relh_adj_tex(0.0);
float water_h_off(0.0), perspective_fovy(0.0), perspective_nclip(0.0), atmosphere(1.0), read_mesh_zmm(0.0);
float light_int_scale(1.0), light_int_scale_l(1.0), snow_depth(0.0), snow_random(0.0);
float init_temperature(DEF_TEMPERATURE);
double camera_zh(0.0);
point mesh_origin(all_zeros), camera_pos(all_zeros);
string user_text;
GLUquadricObj* quadric(0);
lightning l_strike;
colorRGBA bkg_color;
set<unsigned char> keys, keyset;
char game_mode_string[256] = {"640x480:32@85"};
unsigned init_item_counts[] = {2, 2, 2, 6, 6}; // HEALTH, SHIELD, POWERUP, WEAPON, AMMO

// camera variables
double c_radius(DEF_CRADIUS), c_theta(DEF_CTHETA), c_phi(DEF_CPHI), up_theta(DEF_UPTHETA), camera_y(DEF_CAMY);
float sun_rot(0.2), moon_rot(-0.2), light_factor, ball_velocity(15.0), cview_radius(1.0), player_speed(1.0);
vector3d up_vector(0, 1.0, 0);
vector3d cview_dir(all_zeros);
point camera_origin(all_zeros), surface_pos(all_zeros), cpos2;
int orig_window, curr_window;
char player_name[MAX_CHARS] = "Player";


extern bool clear_landscape_vbo, use_dense_voxels;
extern int flight, DISABLE_WATER, DISABLE_SCENERY, camera_invincible, onscreen_display;
extern int tree_coll_level, GLACIATE, UNLIMITED_WEAPONS, destroy_thresh, MAX_RUN_DIST;
extern unsigned NPTS, NRAYS, LOCAL_RAYS, NUM_THREADS, grass_density, max_unique_trees, shadow_map_sz;
extern float fticks, team_damage, self_damage, player_damage, smiley_damage, smiley_speed, tree_deadness, lm_dz_adj, nleaves_scale;
extern float mesh_scale, mesh_scale2, mesh_height_scale, smiley_acc, hmv_scale, last_temp, grass_length, grass_width;
extern point hmv_pos;
extern int coll_id[];
extern vector<bbox> team_starts;
extern player_state *sstates;
extern pt_line_drawer obj_pld;
extern tree_cont_t t_trees;


void init_keyset();
int read_coll_objects(const char *coll_obj_file);
int load_config(string const &config_file);

bool export_modmap(string const &filename);
void reset_planet_defaults();

void verify_wmode(player_state &sstate);

void create_sin_table();



bool check_gl_error(unsigned loc_id) {

	int const error(glGetError());

	if (error) {
		cout << "GL Error " << error << " at location id " << loc_id << ": " << gluErrorString(error) << "." << endl;
		assert(0);
		return 1;
	}
	return 0;
}


void doRedraw(int arg) {

	glutPostWindowRedisplay(curr_window);
}


void clear_context() {

	reset_textures();
	free_dlists();
	clear_shaders();
	reset_snow_vbos();
	update_grass_vbos();
	clear_tree_vbos();
	reset_tiled_terrain_state();
	clear_landscape_vbo = 1;
}


void init_context() {

	screen_reset = 1;
	glFinish();
}


void maximize() {

	//glutHideWindow();
	clear_context();
	glutGameModeString(game_mode_string);
	curr_window   = glutEnterGameMode();
	init_window();
	glFinish();
	ww2           = window_width;
	wh2           = window_height;
	window_width  = gmww;
	window_height = gmwh;
	init_context();
	glutSetCursor(GLUT_CURSOR_NONE);
	nop_frame     = 1;
}


void un_maximize() {

	//glutShowWindow();
	clear_context();
	glutLeaveGameMode();
	curr_window   = orig_window;
	window_width  = ww2;
	window_height = wh2;
	init_context();
	//nop_frame = 1;
}


void init_window() {

	glutSetCursor(GLUT_CURSOR_CROSSHAIR);
	glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard2);

	if (KBD_HANDLER) {
		glutIgnoreKeyRepeat(1);
		glutKeyboardUpFunc(keyboard_up);
		glutSpecialUpFunc(keyboard2_up);
		init_keyset();
	}
	glutPassiveMotionFunc(mousePassiveMotion);

    // Initialize GL
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	glClearColor(0.0, 0.0, 0.0, 0.0);
}


void enable_blend() {

	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);
}


void disable_blend() {

	glDisable(GL_BLEND);
	glDisable(GL_POINT_SMOOTH);
	//glDisable(GL_LINE_SMOOTH);
	//glDisable(GL_POLYGON_SMOOTH);
}


void set_array_client_state(bool va, bool tca, bool na, bool ca) {

	bool const enables[4] = {va, tca, na, ca};
	int  const arrays [4] = {GL_VERTEX_ARRAY, GL_TEXTURE_COORD_ARRAY, GL_NORMAL_ARRAY, GL_COLOR_ARRAY};

	for (unsigned i = 0; i < 4; ++i) {
		if (enables[i]) {
			glEnableClientState(arrays[i]);
		}
		else {
			glDisableClientState(arrays[i]);
		}
	}
}


void set_lighted_sides(int num) {

	assert(num == 1 || num == 2);
	float lmodel_side[1];
	lmodel_side[0] = float(num - 1);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_side);
}


void enable_point_specular() {

	float lmodel_localviewer[] = {1};
	glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_localviewer);
}


void disable_point_specular() {

	float lmodel_localviewer[] = {0};
	glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_localviewer);
}


void set_light_atten(int light, float attenuation) {

	float vals[4];
	glGetLightfv(light, GL_POSITION, vals);

	if (vals[3] != 0.0) { // point light
		glLightf(light, GL_CONSTANT_ATTENUATION,  attenuation);
		glLightf(light, GL_LINEAR_ATTENUATION,    0.0);
		glLightf(light, GL_QUADRATIC_ATTENUATION, 0.0);
	}
	else { // directional light
		int params[2] = {GL_AMBIENT, GL_DIFFUSE};

		for (unsigned i = 0; i < 2; ++i) {
			glGetLightfv(light, params[i], vals);
			for (unsigned d = 0; d < 4; ++d) vals[d] /= attenuation;
			glLightfv(light, params[i], vals);
		}
	}
}


void reset_fog() {

	glFogi(GL_FOG_MODE, GL_LINEAR);
	glFogfv(GL_FOG_COLOR, (float *)&GRAY);
	glFogf(GL_FOG_DENSITY, 0.2);
	glFogf(GL_FOG_START, 0.1);
	glFogf(GL_FOG_END, ((world_mode == WMODE_INF_TERRAIN) ? get_inf_terrain_fog_dist() : 2.5*Z_SCENE_SIZE));
	glFogf(GL_FOG_INDEX, 0.0);
}


void set_gl_params() {

	float lmodel_ambient[]     = {0.0, 0.0, 0.0, 1.0};
    float lmodel_localviewer[] = {0.0};
	float lmodel_side[]        = {0.0};
	setup_fog_coord();
	reset_fog();
	glDepthFunc(GL_LESS);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_localviewer);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_side);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glDisable(GL_COLOR_MATERIAL);
	set_specular(0.0, 1.0);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);

	glHint(GL_FOG_HINT, GL_NICEST); // doesn't work correctly on my laptop
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); // GL_FASTEST
	glHint(GL_LINE_SMOOTH_HINT,    GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT,   GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
}


void reset_camera_pos() {

	if (world_mode == WMODE_UNIVERSE) return;
	camera_origin = mesh_origin;
	up_theta      = DEF_UPTHETA;
	c_radius      = DEF_CRADIUS;
	c_theta       = DEF_CTHETA;
	c_phi         = DEF_CPHI;
	camera_y      = DEF_CAMY;
	if (world_mode == WMODE_UNIVERSE && !island) camera_origin.z = zcenter;
	if (passive_motion) m_button = GLUT_LEFT_BUTTON;
}


void set_perspective(float fovy, float nc_scale) {

	perspective_fovy  = fovy;
	perspective_nclip = nc_scale*NEAR_CLIP;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy, ((GLdouble)window_width)/window_height, perspective_nclip, FAR_CLIP);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


void check_zoom() {

	float fovy(PERSP_ANGLE);

	if (do_zoom == 1) {
		fovy /= ZOOM_FACTOR;
		do_zoom = 2;
	}
	else if (do_zoom == 2) {
		do_zoom = 1;
	}
	set_perspective(fovy, 1.0);
}


void check_xy_offsets() {

	if (camera_view) return;

	while (xoff >= MAX_RUN_DIST) {
		xoff          -= MAX_RUN_DIST;
		surface_pos.x -= MAX_RUN_DIST*DX_VAL;
	}
	while (xoff <= -MAX_RUN_DIST) {
		xoff          += MAX_RUN_DIST;
		surface_pos.x += MAX_RUN_DIST*DX_VAL;
	}
	while (yoff >= MAX_RUN_DIST) {
		yoff          -= MAX_RUN_DIST;
		surface_pos.y -= MAX_RUN_DIST*DY_VAL;
	}
	while (yoff <= -MAX_RUN_DIST) {
		yoff          += MAX_RUN_DIST;
		surface_pos.y += MAX_RUN_DIST*DY_VAL;
	}
}


float calc_speed() {

	float speed(pow(3.0f, min(2, do_run)));
	if (camera_in_air) speed *= CAMERA_AIR_CONT; // what about smileys?
	if (!KBD_HANDLER)  speed /= 0.75;
	return speed;
}


void update_cpos() {

	if (world_mode != WMODE_UNIVERSE) {
		if (fabs(c_phi) < 0.001 || fabs(c_phi - PI) < 0.001 || fabs(c_phi - TWO_PI) < 0.001) c_phi += 0.01;
		
		if (!spectate) {
			cview_dir = rtp_to_xyz(1.0, c_theta, c_phi); // spherical coordinates
			cview_dir.negate();
		}
		cview_radius = c_radius;
	}
	get_camera_pos_v2(camera_pos);
}


void move_camera_pos(vector3d const &v, float dist) { // remember that dist is negative

	if (!camera_surf_collide || flight) {
		surface_pos += v*dist;
	}
	else { // normal ground movement - should speed depend on orientation or not?
		float const xy_scale(dist*(v.mag()/v.xy_mag()));
		surface_pos.x += xy_scale*v.x;
		surface_pos.y += xy_scale*v.y;
	}
}


void advance_camera(int dir) {

	advanced = 1;

	if (world_mode == WMODE_UNIVERSE) { // universe
		bool const hyperspeed(do_run == 2);
		if (player_ship_inited()) player_ship().thrust(dir, speed_mult, hyperspeed);
	}
	else if (world_mode != WMODE_UNIVERSE) {
		if (camera_mode != 1 || (map_mode && world_mode != WMODE_INF_TERRAIN)) return;
		vector3d v;
		float dist(fticks*speed_mult*player_speed*GROUND_SPEED*calc_speed());
		if (game_mode && sstates != NULL) dist *= sstates[CAMERA_ID].get_rspeed_scale();
		
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
}


void change_terrain_zoom(float val) {

	if (island) {
		mesh_scale     = 1.0;
		mesh_scale2   /= val;
		rand_gen_index = rand();
		regen_trees(1, 1);
		compute_volume_matrix(); // make lightning strike the new tree(s)
		return;
	}
	last_temp         = -100.0; // force update
	mesh_scale2       = 1.0;
	camera_change     = 1;
	mesh_scale_change = 1;
	update_mesh(val, 2);
	clear_tiled_terrain();
	if (world_mode == WMODE_INF_TERRAIN) regen_trees(1, 1);
	compute_volume_matrix(); // make lightning strike the new tree(s)
	calc_watershed();
}


void change_world_mode() { // switch terrain mode: 0 = normal, 1 = planet, 2 = no redraw, 3 = dynamic terrain

	if (map_mode || universe_only || (disable_universe && disable_inf_terrain)) return;
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
	
	if (world_mode == WMODE_INF_TERRAIN && island) world_mode = WMODE_GROUND; // unsupported
	
	if (world_mode == WMODE_UNIVERSE) { // save camera position parameters
		xoff_ = xoff; yoff_ = yoff; xoff2_ = xoff2; yoff2_ = yoff2;
		camera_pos_ = camera_pos; // fix_player_upv()?
	}
	if (!map_mode) reset_offsets(); // ???
	init_x        = 1;
	star_init     = 0;
	camera_change = 1;
	reset_fog();
	clear_tiled_terrain();
	obj_pld.free_mem();
	//if (world_mode == WMODE_GROUND || world_mode == WMODE_INF_TERRAIN) calc_visibility(SUN_SHADOW | MOON_SHADOW);
	glDrawBuffer(GL_BACK);
	glutPostWindowRedisplay(curr_window);
	
	if (world_mode == WMODE_INF_TERRAIN) {
		mt2        = mesh_type;
		swap(resolution, res_old);
	}
	else if (world_mode == WMODE_GROUND) {
		mesh_type  = mt2;
		swap(resolution, res_old);
		if (combined_gu) regen_trees(1, 0);
	}
}


// This function is called whenever the window is resized. 
// Parameters are the new dimentions of the window
void resize(int x, int y) {

	if (glutGetWindow() != curr_window) return; // only process the current window

	if (init_resize) {
		init_resize = 0;
	}
	else {
		add_uevent_resize(x, y);
	}
    glViewport(0, 0, x, y);
    window_width  = x;
    window_height = y;
	set_perspective(PERSP_ANGLE, 1.0);
	set_gl_params();
	calc_viewing_cone();
	curr_window = glutGetWindow();
	glutPostWindowRedisplay(curr_window);
}


// This function is called whenever the mouse is pressed or released
// button is a number 0 to 2 designating the button
// state is 1 for release 0 for press event
// x and y are the location of the mouse (in window-relative coordinates)
// Note: Function variables are copied to global versions.
void mouseButton(int button, int state, int x, int y) {

	bool const fire_button((button == GLUT_RIGHT_BUTTON || (passive_motion && button == GLUT_LEFT_BUTTON)));
	add_uevent_mbutton(button, state, x, y);

	if ((camera_mode == 1 || world_mode == WMODE_UNIVERSE) && fire_button) {
		b2down = !state;
		return;
	}
	m_button     = button;
	last_mouse_x = x;
	last_mouse_y = y;
	mouse_state  = state;
}


//This function is called whenever the mouse is moved with a mouse button held down.
// x and y are the location of the mouse (in window-relative coordinates)
void mouseMotion(int x, int y) {

	int button(m_button);

	if (screen_reset) {
		last_mouse_x = x;
		last_mouse_y = y;
		screen_reset = 0;
		return;
	}
	add_uevent_mmotion(x, y);
	float dx(float(x - last_mouse_x)), dy(float(y - last_mouse_y));
	if (camera_mode == 1 && passive_motion) button = GLUT_LEFT_BUTTON;

	switch (button) {
	case GLUT_LEFT_BUTTON: // h: longitude, v: latitude
		if (dx == 0 && dy == 0) break;

		if (world_mode == WMODE_UNIVERSE) {
			vector3d delta(dx, dy, 0.0); // mouse delta
			delta *= (1280.0/window_width); // ???
			if (player_ship_inited()) player_ship().turn(delta);
			update_cpos();
			break;
		}
		{
			float c_phi2(c_phi - MOUSE_ANG_ADJ*dy);
			if (camera_mode && world_mode != WMODE_UNIVERSE) c_phi2 = max(0.01f, min((float)PI-0.01f, c_phi2)); // walking on ground
			
			if (dy > 0.0) { // change camera y direction when camera moved through poles and jump over poles (x=0,z=0) to eliminate "singularity"
				if (c_phi2 < 0.0 || (c_phi > PI && c_phi2 < PI)) camera_y *= -1.0;
				if (fabs(c_phi2) < MA_TOLERANCE || fabs(c_phi2 - PI) < MA_TOLERANCE) dy += 1.0;
			}
			else if (dy < 0.0) {
				if (c_phi2 > TWO_PI || (c_phi < PI && c_phi2 > PI)) camera_y *= -1.0;
				if (fabs(c_phi2 - TWO_PI) < MA_TOLERANCE || fabs(c_phi2 - PI) < MA_TOLERANCE) dy -= 1.0;
			}
			c_theta -= MOUSE_ANG_ADJ*dx*camera_y;
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
		break;

	case GLUT_RIGHT_BUTTON: // v: radius, h: up_vector
		if (camera_mode == 1) break;
		//if (display_mode & 0x80) {mouse_draw_on_ground(x, y);break;}

		if (!camera_view) {
			up_theta += MOUSE_ANG_ADJ*dx;
			up_theta  = fix_angle(up_theta);
			c_radius  = c_radius*(1.0 + MOUSE_R_ADJ*dy);
			if (c_radius <= 0.05*MOUSE_R_ADJ) c_radius = 0.05*MOUSE_R_ADJ;
			update_cpos();
		}
		break;
	}
	last_mouse_x = x;
	last_mouse_y = y;

	if (passive_motion) { // wrap the pointer when it goes off screen
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
	if (dx != 0.0 || dy != 0.0) glutPostWindowRedisplay(curr_window); // Schedule a new display event
}


void mousePassiveMotion(int x, int y) {

	if (passive_motion) mouseMotion(x, y);
}


// This function is called whenever there is a keyboard input
// key is the ASCII value of the key pressed (esc = 27, enter = 13, backspace = 8, tab = 9, del = 127)
// x and y are the location of the mouse
void keyboard_proc(unsigned char key, int x, int y) {

	int mtime2;
	static int lmtype(0);

    switch (key) {
	case 0x1B: // ESC key (27)
		cout << "quitting" << endl;
		exit(0); // quit
		break;
	
	case 'm': // maximize/minimize
		if (!displayed) break;
		mtime2 = glutGet(GLUT_ELAPSED_TIME);
		if (min_time != 0 && (mtime2 - min_time) < MIN_TIME_MS) break;
		min_time = mtime2;
		if (maximized) un_maximize(); else maximize();
		maximized = !maximized;
		displayed = 0;
		break;

	case 'r': // run mode (always move forward)
		run_forward = !run_forward;
		break;
	case 'V': // change mouse mode
		passive_motion = !passive_motion;
		break;

	case 'C': // recreate mesh / add red ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_RED); // red
			break;
		}
		if (mesh_seed != 0 || read_heightmap) break;
		if (mesh_type != lmtype) {
			star_init = 0;
			lmtype = mesh_type;
		}
		rand_gen_index = rand();
		gen_scene(1, (world_mode == WMODE_GROUND), 0, 0, 0);
		regen_lightmap();
		recreated     = 1;
		camera_change = 1;
		break;

	case 'x': // toggle animation
		animate = !animate;
		if (animate) reset_timing = 1;
		break;
	case 't': // animation - movement (freeze frame on objects), show star streams in universe mode
		animate2 = !animate2;
		if (animate2) reset_timing = 1;
		break;

	case 'b': // begin motion animation
		begin_motion = !begin_motion;
		free_dodgeballs(1, 1);
		break;

	case 'y': // save eventlist
		save_ueventlist();
		break;

	case 'k': // change drawing model for mesh (filled polygons vs. wireframe)
		if (map_mode) {
			map_color = !map_color;
		}
		else {
			draw_model = !draw_model;
			set_fill_mode();
		}
		break;

	case 'i': // change mesh type (flat/mountain island vs. mountains) / toggle autopilot
		if (world_mode == WMODE_UNIVERSE) {
			toggle_autopilot();
			break;
		}
		mesh_type  = (mesh_type+1)%3;
		mesh_scale = mesh_scale2 = 1.0;
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
		if (map_mode) {
			map_zoom /= MAP_ZOOM;
		}
		else if (world_mode == WMODE_UNIVERSE) {
			onscreen_display = !onscreen_display;
		}
		else {
			do_zoom = !do_zoom;
			calc_viewing_cone();
		}
		break;

	case 'f': // print framerate and stats
		show_framerate = 1;
		break;
	case 'g': // pause/resume playback of eventlist
		pause_frame = !pause_frame;
		break;
	case 'G': // toggle show framerate (used to be z)
		display_framerate = !display_framerate;
		break;

	case 'p': // reset camera / change fire primary
		if (world_mode == WMODE_UNIVERSE) {
			change_fire_primary();
			break;
		}
		reset_camera_pos();
		break;

	case 'v': // reset camera and change camera mode from air to surface
		reset_camera_pos();
		if (world_mode != WMODE_GROUND) break; // universe mode
		gamemode_rand_appear();
		camera_mode   = !camera_mode;
		camera_reset  = 1;
		camera_change = 1;
		if (camera_mode == 1) camera_invincible = 1; // in air (else on ground)
		break;

	case 'h': // change camera surface collision detection
		if (world_mode == WMODE_UNIVERSE) {
			univ_planet_lod = !univ_planet_lod;
			break;
		}
		camera_surf_collide = !camera_surf_collide;
		camera_change       = 1;
		break;

	case 'j': // smooth camera collision detection / hold fighters
		if (world_mode == WMODE_UNIVERSE) {
			toggle_hold_fighters();
			break;
		}
		camera_coll_smooth = !camera_coll_smooth;
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
		if (world_mode == WMODE_UNIVERSE) player_ship().switch_weapon(1); else switch_player_weapon(-1);
		break;
	case 'e': // next weapon
		if (world_mode == WMODE_UNIVERSE) player_ship().switch_weapon(0); else switch_player_weapon(1);
		break;

		// shadows
	case 'o':
		shadow_detail = min(6, shadow_detail+1);
		cout << "shadow detail = " << shadow_detail << endl;
		break;
	case 'u':
		shadow_detail = max(0, shadow_detail-1);
		cout << "shadow detail = " << shadow_detail << endl;
		break;

	case '=': // increase temp
		temperature += TEMP_INCREMENT;
		cout << "Temperature = " << temperature << " degrees C." << endl;
		temp_change = 1;
		break;
	case '-': // decrease temp
		temperature -= TEMP_INCREMENT;
		temperature  = max(temperature, (island ? 1.0f : ABSOLUTE_ZERO));
		cout << "Temperature = " << temperature << " degrees C." << endl;
		temp_change = 1;
		break;

	case '\'': // increase timestep
		change_timestep(D_TIMESTEP);     break;
	case ';': // decrease timestep
		change_timestep(1.0/D_TIMESTEP); break;

	// resolution
	case ',':
		++resolution;
		break;
	case '.':
		if (resolution > 1) --resolution;
		break;

	case 'T': // delete and generate tree(s) and scenery
		if (world_mode == WMODE_UNIVERSE) {
			player_ship().reset_ammo();
		}
		else {
			rand_gen_index = rand();
			//gen_scenery();
			//regen_trees(1, 0);
			//compute_volume_matrix(); // make lightning strike the new tree(s)
			gen_scene(0, 1, 1, 0, 1);
		}
		break;

	case 'R': // run mode
		++do_run;
		if (do_run > 3) do_run = 0;
		if (world_mode == WMODE_UNIVERSE) change_speed_mode(do_run);
		cout << "run mode = " << do_run << endl;
		break;

	case 'K': // toggle overhead map mode
		if (map_mode) map_mode = 0; else map_mode = 2;
		break;

	case 'U':
		if (world_mode != WMODE_UNIVERSE) { // toggle universe background mode
			combined_gu = !combined_gu;
			if (combined_gu) setup_current_system(); else reset_planet_defaults(); // have to do this so that regen_trees gets correct vegetation
			regen_trees(1, 0);
			reset_tiled_terrain_state();
			
			if (!combined_gu) {
				setup_landscape_tex_colors(ALPHA0, ALPHA0);
				glDisable(get_universe_ambient_light());
				calc_visibility(SUN_SHADOW); // reclaculate sun
				DISABLE_WATER = INIT_DISABLE_WATER;
			}
			create_landscape_texture();
		}
		break;

	// leaf colors
	case 'Q':
		leaf_color_coherence -= 0.05;
		leaf_color_changed    = 1;
		break;
	case 'W':
		leaf_color_coherence += 0.05;
		leaf_color_changed    = 1;
		break;
	case 'A':
		leaf_base_color.red  -= 0.1;
		leaf_color_changed    = 1;
		break;
	case 'S':
		if (world_mode == WMODE_UNIVERSE) {toggle_player_ship_stop(); break;}
		leaf_base_color.red  += 0.1;
		leaf_color_changed    = 1;
		break;
	case 'Z':
		if (map_mode) {
			map_zoom *= MAP_ZOOM;
		}
		else {
			leaf_base_color.green -= 0.1;
		}
		break;
	case 'X':
		leaf_base_color.green += 0.1; break;
	case 'O':
		tree_color_coherence  -= 0.1; break;
	case 'P':
		tree_color_coherence  += 0.1; break;
	case 'E':
		leaf_color_coherence  = 0.5;
		tree_color_coherence  = 0.2;
		leaf_base_color.red   = 0.2;
		leaf_base_color.green = 1.0;
		leaf_color_changed    = 1;
		break;

	case 'L': // increase terrain zoom
		if (mesh_seed != 0 || read_heightmap) break; // ???
		if (!world_mode == WMODE_UNIVERSE) change_terrain_zoom(2.0);
		break;
	case 'Y': // decrease terrain zoom
		if (mesh_seed != 0 || read_heightmap) break; // ???
		if (!world_mode == WMODE_UNIVERSE) change_terrain_zoom(0.5);
		break;

	// object enables
	case 'c': // toggle smileys / add blue ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_BLUE); // blue
			break;
		}
		free_dodgeballs(0, 1);
		obj_groups[coll_id[SMILEY]].toggle_enable();
		if (obj_groups[coll_id[SMILEY]].enabled) init_smileys();
		break;
	case 'B': // precipitation / add neutral ships
		if (world_mode == WMODE_UNIVERSE) {
			add_other_ships(ALIGN_NEUTRAL); // neutral
			break;
		}
		obj_groups[coll_id[PRECIP]].toggle_enable();
		break;

	case 'N': // decrease precipitation rate by 1.5X
		obj_groups[coll_id[PRECIP]].update_app_rate(1.0/1.5, 2, 1000);
		obj_pld.free_mem();
		break;
	case 'M': // increase precipitation rate by 1.5X
		obj_groups[coll_id[PRECIP]].update_app_rate(1.5, 2, 1000);
		break;

	case 'H': // save mesh state/modmap
		if (world_mode != WMODE_UNIVERSE) save_state(state_file); else export_modmap("output.modmap");
		break;
	case 'J': // load mesh state
		if (world_mode != WMODE_UNIVERSE) load_state(state_file);
		break;
	case 'I': // write mesh points
		if (world_mode != WMODE_UNIVERSE) write_mesh(mesh_file);
		break;

	// screenshots
	case 'D': // .raw
		cout << "screenshot.raw: width = " << window_width << ", height = " << window_height << endl;
		glReadBuffer(GL_FRONT);
		screenshot(window_width, window_height, "./");
		break;

		// *** WHAT ABOUT PNG? ***
	case 'F': // .jpg
		cout << "screenshot.jpg: width = " << window_width << ", height = " << window_height << endl;
		glReadBuffer(GL_FRONT);
		write_jpeg(window_width, window_height, "./");
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

	case ' ': // fire key
		fire_weapon();
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

	case '1': // toggle mesh draw
		display_mode ^= 0x01;   break;
	case '2': // toggle grass/snow draw
		display_mode ^= 0x02;   break;
	case '3': // toggle water/ice
		display_mode ^= 0x04;   break;
	case '4': // toggle occlusion culling
		display_mode ^= 0x08;   break;
	case '5': // toggle extended/tiled mesh draw (was debugging)
		display_mode ^= 0x10;   break;
	case '6': // toggle water reflections
		display_mode ^= 0x20;   break;
	case '7': // toggle snow accumulation and clouds
		display_mode ^= 0x40;   break;
	case '8': // toggle water caustics/smoke (change text draw mode - stroke vs. bitmap - disabled)
		display_mode ^= 0x80;   break;
	case '9': // toggle ocean waves and leaf wind
		display_mode ^= 0x0100; break;
	case '0': // change shadow model (stencil shadows and overhead map shadows)
		if (world_mode == WMODE_UNIVERSE) {
			univ_stencil_shadows = !univ_stencil_shadows;
		}
		else {
			use_stencil_shadows  = !use_stencil_shadows;
		}
		break;

	case '\\': // enable dynamic particles (to test dynamic lighting, dynamic shadows, and collision detection)
		display_mode ^= 0x0200;
		d_part_sys.clear();
		break;
	}
	glutPostWindowRedisplay(curr_window); // Schedule a new display event
}


void print_wind() {

	cout << "wind: "; wind.print(); cout << endl;
}


void keyboard2(int key, int x, int y) {

	add_uevent_keyboard_special(key, x, y);

	switch (key) {
	case GLUT_KEY_UP:
		if (map_mode) {
			map_y += int(map_zoom*MAP_SHIFT) + 1;
		}
		else {
			wind.y += WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_DOWN:
		if (map_mode) {
			map_y -= int(map_zoom*MAP_SHIFT) + 1;
		}
		else {
			wind.y -= WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_LEFT:
		if (map_mode) {
			map_x -= int(map_zoom*MAP_SHIFT) + 1;
		}
		else {
			wind.x -= WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_RIGHT:
		if (map_mode) {
			map_x += int(map_zoom*MAP_SHIFT) + 1;
		}
		else {
			wind.x += WIND_ADJUST;
			print_wind();
		}
		break;

	case GLUT_KEY_F1: // switch terrain mode: 0 = normal, 1 = planet, 2 = no redraw, 3 = dynamic terrain
		change_world_mode();
		break;

	case GLUT_KEY_F2: // switch game mode
		if (world_mode == WMODE_UNIVERSE) break;
		game_mode = (game_mode + 1)%3;
		change_game_mode();
		break;

	case GLUT_KEY_F4: // switch weapon mode
		if (sstates != NULL) {
			++sstates[CAMERA_ID].wmode;
			sstates[CAMERA_ID].verify_wmode();
		}
		break;

	case GLUT_KEY_F5: // toggle large/small trees
		if (world_mode != WMODE_GROUND || (num_trees == 0 && t_trees.empty())) break;
		tree_mode = (tree_mode+1)%4; // 0=none, 1=large, 2=small, 3=large+small
		gen_scene(0, 1, 1, 0, 1);
		break;

	case GLUT_KEY_F6: // toggle draw_mesh shader
		draw_mesh_shader ^= 1;
		break;

	case GLUT_KEY_F7: // toggle auto time advance
		obj_groups[coll_id[PRECIP]].app_rate = 40;
		auto_time_adv = (auto_time_adv+1)%5;
		cout << "Auto time advance = " << auto_time_adv << endl;
		break;

	case GLUT_KEY_F8: // toggle spectator gameplay mode
		if (!spectate && (num_smileys == 0 || !obj_groups[coll_id[SMILEY]].enabled)) break;
		if (spectate) {camera_reset = camera_change = 1;}
		spectate = !spectate;
		break;

	case GLUT_KEY_F9: // switch to fullscreen mode
		glutFullScreen();
		break;
	}
	glutPostWindowRedisplay(curr_window); // Schedule a new display event
}


void init_keyset() {

	string keyvals = "wsad "; // movement

	for (unsigned i = 0; i < keyvals.size(); ++i) {
		keyset.insert(keyvals[i]);
	}
}


unsigned char get_key_other_case(unsigned char key) {

	unsigned char key2(key); // check if shift key was pressed while holding down a key
	if (key >= 'a' && key <= 'z') key2 = key + ('A' - 'a');
	if (key >= 'A' && key <= 'Z') key2 = key + ('a' - 'A');
	return key2;
}


void keyboard_up(unsigned char key, int x, int y) {

	if (!KBD_HANDLER || kbd_text_mode || key == 13) return; // ignore text mode and enter keys
	if (keyset.find(key) != keyset.end()) add_uevent_keyboard_up(key, x, y);
	keyset_it it(keys.find(key));

	if (it == keys.end()) {
		unsigned char key2(get_key_other_case(key));
		if (key2 != key) it = keys.find(key2);

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

	// nothing
}


void exec_text(string const &text) {

	if (world_mode == WMODE_UNIVERSE) { // handled in universe_control.cpp
		exec_universe_text(text);
		return;
	}
	cout << "Text: " << text << endl;
	print_text_onscreen(text.c_str(), WHITE, 1.0, TICKS_PER_SECOND, 999);
}


void keyboard(unsigned char key, int x, int y) {

	add_uevent_keyboard(key, x, y);

	if (key == 13) { // enter key - toggle text mode
		if (kbd_text_mode) exec_text(user_text);
		user_text.clear();
		kbd_text_mode = !kbd_text_mode;
		return;
	}
	if (kbd_text_mode) {
		if (key == 8) { // delete key
			if (!user_text.empty()) user_text.erase(user_text.begin()+user_text.size()-1); // pop_back() for string
			return;
		}
		user_text.push_back(key);
		return;
	}
	if (!KBD_HANDLER) {
		keyboard_proc(key, x, y);
		return;
	}
	if (keys.find(key) != keys.end()) { // can happen with control clicks
		cout << "Warning: Keyboard event for key " << key << " (" << int(key) << ") which has alredy been pressed." << endl;
		//assert(0); // too strong?
		return;
	}
	keys.insert(key);
	if (keyset.find(key) == keyset.end()) keyboard_proc(key, x, y);
}


void proc_kbd_events() {

	if (!KBD_HANDLER) return;

	for (keyset_it it = keys.begin(); it != keys.end(); ++it) {
		if (keyset.find(*it) != keyset.end()) keyboard_proc(*it, 0, 0); // x and y = ?
	}
}


inline bool read_int  (FILE *fp, int      &val) {return (fscanf(fp, "%i", &val) == 1);}
inline bool read_uint (FILE *fp, unsigned &val) {return (fscanf(fp, "%u", &val) == 1);}
inline bool read_float(FILE *fp, float    &val) {return (fscanf(fp, "%f", &val) == 1);}
inline bool read_str  (FILE *fp, char     *val) {return (fscanf(fp, "%s",  val) == 1);}

inline bool read_bool (FILE *fp, bool     &val) {
	int tmp;
	if (fscanf(fp, "%i", &tmp) != 1) return 0;
	val = (tmp != 0);
	return 1;
}


int check_for_config_defaults(const char *def_file) {

	assert(def_file != NULL);
	string config_file;
	ifstream in(def_file);
	if (!in.good()) return 0;

	while (in >> config_file) {
		if (!config_file.empty() && config_file[0] != '#') { // not commented out
			cout << "Using config file " << config_file << "." << endl;
			load_config(config_file);
		}
	}
	return 1;
}


void fire_weapon() {

	fire_key = 1;

	if (world_mode != WMODE_UNIVERSE) {
		assert(sstates != NULL);
		sstates[CAMERA_ID].gamemode_fire_weapon();
	}
}


bool has_extension(string const &ext) { // is this always correct?

	return (strstr((const char *)glGetString(GL_EXTENSIONS), ext.c_str()) != NULL);
}


bool open_file(FILE *&fp, char const *const fn, string const &file_type, char const *const mode) {

	fp = fopen(fn, mode);
	if (fp != NULL) return 1;
	cout << "Could not open " << file_type << " file '" << fn << "'." << endl;
	return 0;
}


void alloc_if_req(char *&fn, const char *def_fn) {

	if (fn == def_fn) fn = new char[MAX_CHARS];
}


void cfg_err(string const &str, int &error) {

	cout << "Error reading " << str << " from config file." << endl;
	error = 1;
}


// should be moved to another file eventually...
// should use a hashtable here
int load_config(string const &config_file) {

	int gms_set(0), color_bit_depth(32), refresh_rate(75), error(0);
	char strc[MAX_CHARS] = {0}, md_fname[MAX_CHARS] = {0}, we_fname[MAX_CHARS] = {0};
	gmww          = 640;
	gmwh          = 480;
	mesh_draw     = NULL;
	water_enabled = NULL;
	FILE *fp;
	if (!open_file(fp, config_file.c_str(), "input configuration file")) return 0;

	while (read_str(fp, strc)) { // slow but should be OK
		string const str(strc);

		if (str[0] == '#') { // comment
			char letter(getc(fp));
			while (letter != '\n' && letter != EOF && letter != 0) letter = getc(fp);
		}
		else if (str == "grass_density") {
			if (!read_uint(fp, grass_density)) cfg_err("grass density", error);
		}
		else if (str == "grass_size") {
			if (!read_float(fp, grass_length) || !read_float(fp, grass_width) || grass_length <= 0.0 || grass_width <= 0.0) {
				cfg_err("grass size", error);
			}
		}
		else if (str == "ntrees") {
			if (!read_int(fp, num_trees) || num_trees < 0) cfg_err("number of trees", error);
		}
		else if (str == "max_unique_trees") {
			if (!read_uint(fp, max_unique_trees)) cfg_err("max_unique_trees", error);
		}
		else if (str == "nleaves_scale") {
			if (!read_float(fp, nleaves_scale) || nleaves_scale <= 0.0) cfg_err("nleaves_scale", error);
		}
		else if (str == "nsmileys") {
			if (!read_int(fp, num_smileys) || num_smileys < 0) cfg_err("number of smileys", error);
		}
		else if (str == "num_items") {
			for (unsigned n = 0; n < sizeof(init_item_counts)/sizeof(unsigned); ++n) {
				if (!read_uint(fp, init_item_counts[n])) {
					cfg_err("number of items", error);
					break;
				}
			}
		}
		else if (str == "verbose") {
			if (!read_int(fp, verbose_mode)) cfg_err("verbose mode", error);
		}
		else if (str == "game_mode_string") {
			if (fscanf(fp, "%s%i%i", game_mode_string, &gmww, &gmwh) != 3 || gmww <= 0 || gmwh <= 0) cfg_err("game mode string", error);
			gms_set = 1;
		}
		else if (str == "window_width") {
			if (!read_int(fp, gmww) || gmww < 1) cfg_err("gmww command", error);
		}
		else if (str == "window_height") {
			if (!read_int(fp, gmwh) || gmwh < 1) cfg_err("gmwh command", error);
		}
		else if (str == "mesh_size") {
			if (fscanf(fp, "%i%i%i", &MESH_X_SIZE, &MESH_Y_SIZE, &MESH_Z_SIZE) != 3) cfg_err("mesh size command", error);
		}
		else if (str == "scene_size") {
			if (fscanf(fp, "%f%f%f", &X_SCENE_SIZE, &Y_SCENE_SIZE, &Z_SCENE_SIZE) != 3) cfg_err("scene size command", error);
		}
		else if (str == "shadow_map_sz") { // 0 disables
			if (!read_uint(fp, shadow_map_sz)) cfg_err("shadow_map_sz command", error);
		}
		else if (str == "color_bit_depth") {
			if (!read_int(fp, color_bit_depth)) cfg_err("color_bit_depth command", error);
			if (color_bit_depth < 1 || (color_bit_depth >= 8 && color_bit_depth%8 != 0)) {cout << "Invalid color_bit_depth command in config file." << endl; error = 1;}
		}
		else if (str == "refresh_rate") {
			if (!read_int(fp, refresh_rate) || refresh_rate < 1) cfg_err("refresh_rate command", error);
		}
		else if (str == "load_hmv") {
			if (!read_int(fp, load_hmv)) cfg_err("load_hmv command", error);
			if (fscanf(fp, "%f%f%f%f", &hmv_pos.x, &hmv_pos.y, &hmv_pos.z, &hmv_scale) != 4) cfg_err("load_hmv command", error);
		}
		else if (str == "load_coll_objs") {
			if (!read_int(fp, load_coll_objs)) cfg_err("load coll objs command", error);
		}
		else if (str == "glaciate") {
			if (!read_int(fp, GLACIATE)) cfg_err("glaciate command", error);
		}
		else if (str == "dynamic_mesh_scroll") {
			if (!read_int(fp, dynamic_mesh_scroll)) cfg_err("dynamic_mesh_scroll command", error);
		}
		else if (str == "use_waypoints") {
			if (!read_bool(fp, use_waypoints)) cfg_err("use_waypoints command", error);
		}
		else if (str == "group_back_face_cull") {
			if (!read_bool(fp, group_back_face_cull)) cfg_err("group_back_face_cull command", error);
		}
		else if (str == "gravity") {
			if (!read_float(fp, base_gravity)) cfg_err("gravity command", error);
		}
		else if (str == "temperature") {
			if (!read_float(fp, init_temperature)) cfg_err("temperature command", error);
			temperature = init_temperature;
		}
		else if (str == "mesh_seed") {
			if (!read_int(fp, mesh_seed)) cfg_err("mesh seed scale command", error);
		}
		else if (str == "mesh_height") {
			if (!read_float(fp, mesh_height_scale)) cfg_err("mesh height scale command", error);
		}
		else if (str == "mesh_scale") {
			if (!read_float(fp, mesh_scale)) cfg_err("mesh scale command", error);
		}
		else if (str == "mesh_z_cutoff") {
			if (!read_float(fp, mesh_z_cutoff)) cfg_err("mesh z cutoff command", error);
		}
		else if (str == "disabled_mesh_z") {
			if (!read_float(fp, disabled_mesh_z)) cfg_err("disabled mesh z command", error);
		}
		else if (str == "relh_adj_tex") {
			if (!read_float(fp, relh_adj_tex)) cfg_err("relh_adj_tex command", error);
		}
		else if (str == "water_h_off") {
			if (!read_float(fp, water_h_off)) cfg_err("water_h_off command", error);
		}
		else if (str == "lm_dz_adj") {
			if (!read_float(fp, lm_dz_adj) || lm_dz_adj < 0.0) cfg_err("lm_dz_adj command", error);
		}
		else if (str == "wind_velocity") {
			if (fscanf(fp, "%f%f%f", &wind.x, &wind.y, &wind.z) != 3) cfg_err("wind_velocity command", error);
		}
		else if (str == "camera_height") {
			if (fscanf(fp, "%lf", &camera_zh) != 1) cfg_err("camera extra height command", error);
		}
		else if (str == "player_start") {
			if (fscanf(fp, "%f%f%f", &surface_pos.x, &surface_pos.y, &surface_pos.z) != 3) cfg_err("player_start command", error);
		}
		else if (str == "tree_size") {
			if (!read_float(fp, tree_size)) cfg_err("tree size command", error);
			mult_leaf_points_by(tree_size);
		}
		else if (str == "tree_deadness") {
			if (!read_float(fp, tree_deadness)) cfg_err("tree deadness command", error);
		}
		else if (str == "leaf_color") {
			if (fscanf(fp, "%f%f%f%f%f", &leaf_base_color.red, &leaf_base_color.green, &leaf_base_color.blue, &leaf_color_coherence, &tree_color_coherence) != 5) {
				cfg_err("leaf color command", error);
			}
		}
		else if (str == "sun_rot") {
			if (!read_float(fp, sun_rot)) cfg_err("sun_rot command", error);
		}
		else if (str == "moon_rot") {
			if (!read_float(fp, moon_rot)) cfg_err("moon_rot command", error);
		}
		else if (str == "universe_only") {
			if (!read_int(fp, universe_only)) cfg_err("universe only command", error);
		}
		else if (str == "disable_universe") {
			if (!read_int(fp, disable_universe)) cfg_err("disable universe command", error);
		}
		else if (str == "disable_inf_terrain") {
			if (!read_int(fp, disable_inf_terrain)) cfg_err("disable inf terrain command", error);
		}
		else if (str == "player_name") {
			if (!read_str(fp, player_name)) cfg_err("player name", error);
		}
		else if (str == "self_damage") {
			if (!read_float(fp, self_damage)) cfg_err("self damage command", error);
		}
		else if (str == "team_damage") {
			if (!read_float(fp, team_damage)) cfg_err("team damage command", error);
		}
		else if (str == "player_damage") {
			if (!read_float(fp, player_damage)) cfg_err("player damage command", error);
		}
		else if (str == "smiley_damage") {
			if (!read_float(fp, smiley_damage)) cfg_err("smiley damage command", error);
		}
		else if (str == "player_speed") {
			if (!read_float(fp, player_speed)) cfg_err("player speed command", error);
		}
		else if (str == "smiley_speed") {
			if (!read_float(fp, smiley_speed)) cfg_err("smiley speed command", error);
		}
		else if (str == "speed_mult") {
			if (!read_float(fp, speed_mult)) cfg_err("speed mult speed command", error);
		}
		else if (str == "left_handed") {
			if (!read_int(fp, left_handed)) cfg_err("left handed command", error);
		}
		else if (str == "unlimited_weapons") {
			if (!read_int(fp, UNLIMITED_WEAPONS)) cfg_err("unlimited weapons command", error);
		}
		else if (str == "smiley_accuracy") {
			if (!read_float(fp, smiley_acc)) cfg_err("smiley accuracy command", error);
		}
		else if (str == "crater_size") {
			if (!read_float(fp, crater_size)) cfg_err("crater size command", error);
		}
		else if (str == "destroy_thresh") {
			if (!read_int(fp, destroy_thresh)) cfg_err("destroy thresh command", error);
		}
		else if (str == "rand_seed") {
			if (!read_int(fp, srand_param)) cfg_err("random seed command", error);
		}
		else if (str == "disable_water") {
			if (!read_int(fp, INIT_DISABLE_WATER)) cfg_err("disable water command", error);
		}
		else if (str == "disable_scenery") {
			if (!read_int(fp, DISABLE_SCENERY)) cfg_err("disable scenery command", error);
		}
		else if (str == "disable_shaders") {
			if (!read_bool(fp, disable_shaders)) cfg_err("disable shaders command", error);
		}
		else if (str == "use_dense_voxels") {
			if (!read_bool(fp, use_dense_voxels)) cfg_err("use_dense_voxels command", error);
		}
		else if (str == "read_landscape") {
			if (!read_int(fp, read_landscape)) cfg_err("read landscape command", error);
		}
		else if (str == "read_heightmap") {
			if (!read_int(fp, read_heightmap)) cfg_err("read heightmap command", error);
		}
		else if (str == "default_ground_tex") {
			if (!read_int(fp, default_ground_tex)) cfg_err("default ground texture command", error);
		}
		else if (str == "ground_effects_level") {
			if (!read_int(fp, ground_effects_level)) cfg_err("ground_effects_level command", error);
		}
		else if (str == "shadow_detail") {
			if (!read_int(fp, shadow_detail)) cfg_err("shadow detail command", error);
		}
		else if (str == "tree_coll_level") {
			if (!read_int(fp, tree_coll_level)) cfg_err("tree coll level command", error);
		}
		else if (str == "free_for_all") {
			if (!read_int(fp, free_for_all)) cfg_err("free for all command", error);
		}
		else if (str == "teams") {
			if (!read_int(fp, teams) || teams <= 0) cfg_err("teams command", error);
		}
		else if (str == "num_dodgeballs") {
			if (!read_int(fp, num_dodgeballs)) cfg_err("num dodgeballs command", error);
			num_dodgeballs = max(num_dodgeballs, 1); // have to have at least 1
		}
		else if (str == "team_start") {
			bbox bb;
			if (fscanf(fp, "%i%f%f%f%f", &bb.index, &bb.x1, &bb.y1, &bb.x2, &bb.y2) != 5) cfg_err("team start command", error);
			if (bb.index < 0 || bb.index >= teams) {cout << "Error: Illegal team specified in team_start command: " << bb.index << " (" << teams << " teams)." << endl; error = 1;}
			else team_starts.push_back(bb);
		}
		else if (str == "coll_obj_file") {
			alloc_if_req(coll_obj_file, dcoll_obj_file);
			if (!read_str(fp, coll_obj_file)) cfg_err("coll_obj_file command", error);
		}
		else if (str == "state_file") {
			alloc_if_req(state_file, dstate_file);
			if (!read_str(fp, state_file)) cfg_err("state_file command", error);
		}
		else if (str == "mesh_file") {
			float rmz(0.0);
			alloc_if_req(mesh_file, dmesh_file);
			if (fscanf(fp, "%s%f%f%i", mesh_file, &mesh_file_scale, &mesh_file_tz, &do_read_mesh) != 4) cfg_err("mesh_file command", error);
			if (read_float(fp, rmz)) read_mesh_zmm = rmz;
		}
		else if (str == "mh_filename_raw") {
			alloc_if_req(mh_filename_raw, dmh_filename_raw);
			if (fscanf(fp, "%s%f%f", mh_filename_raw, &mesh_file_scale, &mesh_file_tz) != 3) cfg_err("mh_filename_raw command", error);
		}
		else if (str == "mh_filename_bmp") {
			alloc_if_req(mh_filename_bmp, dmh_filename_bmp);
			if (fscanf(fp, "%s%f%f", mh_filename_bmp, &mesh_file_scale, &mesh_file_tz) != 3) cfg_err("mh_filename_bmp command", error);
		}
		else if (str == "dem_filename") {
			alloc_if_req(dem_filename, ddem_filename);
			if (!read_str(fp, dem_filename)) cfg_err("dem_filename command", error);
		}
		else if (str == "dem_raw_out") {
			if (ddem_raw_out == ddem_raw_out) dem_raw_out = new char[MAX_CHARS];
			if (!read_str(fp, dem_raw_out)) cfg_err("dem_raw_out command", error);
		}
		else if (str == "ship_def_file") {
			alloc_if_req(ship_def_file, dship_def_file);
			if (!read_str(fp, ship_def_file)) cfg_err("ship_def_file command", error);
		}

		// lighting
		else if (str == "lighting_file") {
			alloc_if_req(lighting_file, NULL);
			int write_mode(0);
			if (fscanf(fp, "%s%i%f", lighting_file, &write_mode, &light_int_scale) != 3) cfg_err("lighting_file command", error);
			(write_mode ? write_light_file : read_light_file) = 1;
		}
		else if (str == "lighting_file_local") {
			alloc_if_req(lighting_file_l, NULL);
			int write_mode(0);
			if (fscanf(fp, "%s%i%f", lighting_file_l, &write_mode, &light_int_scale_l) != 3) cfg_err("lighting_file_local command", error);
			(write_mode ? write_light_file_l : read_light_file_l) = 1;
		}
		else if (str == "num_light_rays") {
			if (fscanf(fp, "%u%u%u", &NPTS, &NRAYS, &LOCAL_RAYS) != 3) cfg_err("num_light_rays command", error);
		}
		else if (str == "num_threads") {
			if (!read_uint(fp, NUM_THREADS) || NUM_THREADS == 0 || NUM_THREADS > 100) cfg_err("num_threads", error);
		}

		// snow
		else if (str == "num_test_snowflakes") { // in millions
			if (!read_uint(fp, num_snowflakes)) cfg_err("num_test_snowflakes", error);
		}
		else if (str == "snow_depth") {
			if (!read_float(fp, snow_depth) || snow_depth < 0.0) cfg_err("snow_depth command", error);
		}
		else if (str == "snow_random") {
			if (!read_float(fp, snow_random)) cfg_err("snow_random command", error);
		}
		else if (str == "snow_file") {
			alloc_if_req(snow_file, NULL);
			int write_mode(0);
			if (fscanf(fp, "%s%i", snow_file, &write_mode) != 2) cfg_err("snow_file command", error);
			(write_mode ? write_snow_file : read_snow_file) = 1;
		}

		else if (str == "scenery_extent") { // in mesh sizes
			if (!read_uint(fp, scenery_extent)) cfg_err("scenery_extent", error);
		}

		else if (str == "mesh_draw_bmp") {
			if (!read_str(fp, md_fname)) cfg_err("mesh_draw_bmp command", error);
		}
		else if (str == "water_enabled_bmp") {
			if (!read_str(fp, we_fname)) cfg_err("water_enabled_bmp command", error);
		}
		else if (str == "end") {
			break;
		}
		else {
			cout << "Unrecognized keyword in input file: " << str << endl;
			error = 1;
		}
		if (error) {cout << "Parse error in config file." << endl; break;}
	}
	fclose(fp);
	
	if (universe_only) {
		world_mode = WMODE_UNIVERSE;
		assert(!disable_universe);
	}
	/*if (read_landscape) {
		DISABLE_SCENERY = 1;
		DISABLE_WATER   = 1;
		num_trees       = 0;
	}*/
	//if (read_heightmap && dynamic_mesh_scroll) cout << "Warning: read_heightmap and dynamic_mesh_scroll are currently incompatible options as the heightmap does not scroll." << endl;
	DISABLE_WATER = INIT_DISABLE_WATER;
	XY_MULT_SIZE  = MESH_X_SIZE*MESH_Y_SIZE; // for bmp_to_chars() allocation
	if (!gms_set) sprintf(game_mode_string, "%ix%i:%i@%i", gmww, gmwh, color_bit_depth, refresh_rate);
	if (strlen(md_fname) > 0) if (!bmp_to_chars(md_fname, mesh_draw))     error = 1;
	if (strlen(we_fname) > 0) if (!bmp_to_chars(we_fname, water_enabled)) error = 1;
	if (error) exit(1);
	return 1;
}


#if 0
#include <omp.h>
 
 void open_mp_test() {
   int th_id, nthreads;
   #pragma omp parallel private(th_id)
   {
     th_id = omp_get_thread_num();
     printf("Hello World from thread %d\n", th_id);
     #pragma omp barrier
     if ( th_id == 0 ) {
       nthreads = omp_get_num_threads();
       printf("There are %d threads\n",nthreads);
     }
   }

   RESET_TIME;
   float val(0);

   //omp_set_num_threads(4);
   #pragma omp parallel for schedule(static,1) 
   for (int i = 0; i < 10000000; ++i) {
	   val += sinf(123*i) + 1.0;
   }
   PRINT_TIME("For Loop");
   cout << "Val: " << val << endl;
   exit(0);
 }
#else
void open_mp_test() {}
#endif



int main(int argc, char** argv) {

	open_mp_test(); // testing
	cout << "Starting 3DWorld."; cout.flush();
	
    // Initialize GLUT
	cout << "."; cout.flush();
    glutInit(&argc, argv);
	cout << "."; cout.flush();
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
	//glutInitDisplayString("rgb double depth>=16 samples>=4");
	cout << "."; cout.flush();
	orig_window = glutCreateWindow("3D World");
	curr_window = orig_window;
	cout << "."; cout.flush();
	init_glew();
	cout << "."; cout.flush();
	init_window();
	cout << ".GL Initialized." << endl;
	//glutFullScreen();
	if (argc == 2) read_ueventlist(argv[1]);
	int rs(1);

	if (srand_param == 1) {
		rs = glutGet(GLUT_ELAPSED_TIME);
	}
	else if (srand_param != 0) {
		rs = srand_param;
	}
	add_uevent_srand(rs);
	uevent_advance_frame();
	--frame_counter;
	l_strike.enabled = -1;
	quadric = gluNewQuadric();
	
	if (quadric == 0) {
		cout << "Failed to create quadric." << endl;
		return 1;
	}
	create_sin_table();
	gen_gauss_rand_arr();
	set_scene_constants();
	check_for_config_defaults(defaults_file);
	load_textures();
	load_flare_textures(); // Sun Flare
	setup_dlists();
	setup_shaders();
	//cout << "Extensions: " << (char*) glGetString(GL_EXTENSIONS) << endl;

	if (!universe_only) { // universe mode should be able to do without these initializations
		init_objects();
		alloc_matrices();
		t_trees.resize(num_trees);
		init_models();
		init_terrain_mesh();
		gen_scene(1, (world_mode == WMODE_GROUND), 0, 0, 0);
		gen_snow_coverage();
		create_object_groups();
		init_game_state();
		build_lightmap(1);
	}
	glutMainLoop(); // Switch to main loop

	// cleanup - never actually gets here, but could put in quit code
	free_textures();

	if (!universe_only) {
		free_models();
		free_scenery();
		delete_matrices();
	}
    return 0;        
}

