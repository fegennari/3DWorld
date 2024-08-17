// 3D World - Full screen postprocessing effects
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "function_registry.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "transform_obj.h"


extern bool water_is_lava, enable_postproc_recolor;
extern unsigned depth_tid, frame_buffer_RGB_tid;
extern int frame_counter, display_mode, show_fog, camera_coll_id, window_width, window_height, animate2;
extern float NEAR_CLIP, FAR_CLIP, fticks, dist_to_fire_sq, water_plane_z, CAMERA_RADIUS;
extern colorRGBA sun_color;

int depth_buffer_frame(0), color_buffer_frame(0);
float cur_explosion_weight(0.0);
colorRGBA vignette_color(ALPHA0);
sphere_t cur_explosion_sphere;

bool player_is_drowning();
float get_player_drunkenness();


void bind_depth_buffer(unsigned tu_id=0) {
	if (frame_counter != depth_buffer_frame) {depth_buffer_to_texture(depth_tid);} // depth texture is not valid for this frame, create it
	depth_buffer_frame = frame_counter;
	assert(depth_tid > 0);
	bind_texture_tu(depth_tid, tu_id);
}
void bind_frame_buffer_RGB(unsigned tu_id=0) {
	if (frame_counter != color_buffer_frame) {frame_buffer_RGB_to_texture(frame_buffer_RGB_tid);} // FB RGB texture is not valid for this frame, create it
	color_buffer_frame = frame_counter;
	assert(frame_buffer_RGB_tid > 0);
	bind_texture_tu(frame_buffer_RGB_tid, tu_id);
}

void draw_ortho_screen_space_triangle() {

	// setup matrices
	fgMatrixMode(FG_PROJECTION);
	fgPushIdentityMatrix();
	fgOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	fgMatrixMode(FG_MODELVIEW);
	fgPushIdentityMatrix();

	enable_blend();
	glDisable(GL_DEPTH_TEST);
	//draw_tquad(1.0, 1.0, 0.0);
	vert_tc_t const verts[3] = {vert_tc_t(point(-1,-1,0), 0,0), vert_tc_t(point(3,-1,0), 2,0), vert_tc_t(point(-1,3,0), 0,2)};
	draw_verts(verts, 3, GL_TRIANGLES); // supposedly a single triangle has better cache performance
	glEnable(GL_DEPTH_TEST);
	disable_blend();

	restore_prev_mvm_pjm_state();
}

void set_xy_step(shader_t &s) {s.add_uniform_vector2d("xy_step", vector2d(1.0/window_width, 1.0/window_height));}

void setup_depth_tex(shader_t &s, int tu_id) {
	set_xy_step(s);
	s.add_uniform_int("depth_tex", tu_id);
	s.add_uniform_float("znear", NEAR_CLIP);
	s.add_uniform_float("zfar",  FAR_CLIP);
}

void fill_screen_white_and_end_shader(shader_t &s) {
	s.set_cur_color(WHITE);
	draw_ortho_screen_space_triangle();
	s.end_shader();
}

// add God rays as a fullscreen shader pass using the depth texture
void add_god_rays() {

	if (world_mode == WMODE_UNIVERSE) return; // not in universe mode
	if (!is_sun_flare_visible()) return; // sun not visible
	bind_depth_buffer();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("god_rays");
	s.begin_shader();
	s.add_uniform_int("depth_tex", 0);
	s.add_uniform_color("sun_color", sun_color);
	s.add_uniform_vector3d("sun_pos", world_space_to_screen_space(get_sun_pos()));
	s.add_uniform_float("aspect_ratio", float(window_width)/float(window_height));
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

class ssao_state_manager_t {
	unsigned width=0, height=0, tid=0, fbo=0;
public:
	void begin() {
		if (width != (unsigned)window_width || height != (unsigned)window_height) { // size change, clear
			clear();
			width  = window_width;
			height = window_height;
		}
		if (!tid) {
			setup_texture(tid, 0, 0, 0);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
		}
		enable_fbo(fbo, tid, 0);
	}
	void end () {disable_fbo();}
	void bind() {bind_2d_texture(tid);}

	void clear() {
		free_texture(tid);
		free_fbo(fbo);
		width = height = 0; // is this needed?
	}
};

ssao_state_manager_t ssao_state_manager;

void add_ssao() {

	// Note: somewhat works, but looks very bad on windows, grass, and clouds
	bool const USE_SSAO_BLUR = 1; // looks better, but slightly slower
	bind_depth_buffer();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("depth_utils.part+screen_space_ao"); // too blocky, doesn't work on transparent objects, sky/clouds, or grass
	//s.set_frag_shader("depth_utils.part+screen_space_ao_v2"); // black checkerboard patterns, needs blurring
	if (USE_SSAO_BLUR) {s.set_prefix("#define WRITE_COLOR", 1);} // FS
	s.begin_shader();
	setup_depth_tex(s, 0);

	if (USE_SSAO_BLUR) {
		// render SSAO weight to an FBO texture
		ssao_state_manager.begin();
		fill_screen_white_and_end_shader(s);
		ssao_state_manager.end();
		// blur SSAO
		ssao_state_manager.bind();
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("screen_space_ao_blur"); // TODO: 2x 1D blurs?
		s.begin_shader();
		s.add_uniform_int("frame_buffer_tex", 0);
		set_xy_step(s);
	}
	fill_screen_white_and_end_shader(s);
}

void add_color_only_effect(string const &frag_shader, float intensity=1.0, float time_scale=1.0, float pos_scale=0.0, colorRGBA const &color_mod=WHITE) {

	static float time(0.0);
	if (animate2) {time += time_scale*fticks;}
	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader(frag_shader);
	s.begin_shader();
	s.add_uniform_float("intensity", intensity);
	if (pos_scale != 0.0) {s.add_uniform_float("pos_scale", pos_scale);} // not all shaders have this uniform
	s.add_uniform_int("frame_buffer_tex", 0);
	s.add_uniform_float("time", time); // may not be used
	s.add_uniform_color("color_mod", color_mod); // may not be used
	select_texture(NOISE_TEX, 1);
	s.add_uniform_int("noise_tex", 1); // Note: used for heat waves effect, could be used for others
	set_xy_step(s); // may not be used
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

void add_vignette(colorRGBA const &color) {

	if (color.A == 0.0) return;
	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("vignette");
	s.begin_shader();
	s.add_uniform_int("frame_buffer_tex", 0);
	s.add_uniform_color("edge_color", color);
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

void postproc_convert_to_grayscale(unsigned xsize, unsigned ysize) {

	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("convert_to_grayscale");
	s.begin_shader();
	// since the screen resolution may be different, we have to scale the texture coordinates
	s.add_uniform_float("xscale", float(xsize)/float(window_width ));
	s.add_uniform_float("yscale", float(ysize)/float(window_height));
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

void add_sphere_refract_effect(sphere_t const &sphere, float intensity) {

	if (intensity == 0.0) return;
	//static float time(0.0); if (animate2) {time += fticks;}
	point const center(world_space_to_screen_space(sphere.pos));
	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("sphere_refract_screen");
	s.set_prefix("#define CHROMATIC_REFRACT", 1); // FS
	s.begin_shader();
	s.add_uniform_int("frame_buffer_tex", 0);
	//s.add_uniform_float("time", time); // unused in shader
	s.add_uniform_float("aspect_ratio", float(window_width)/float(window_height));
	s.add_uniform_float("intensity", CLIP_TO_01(intensity)); // 1.0 at T=0, 0.0 at T=1
	s.add_uniform_float("radius", sphere.radius/p2p_dist(get_camera_pos(), sphere.pos)); // divide distance/depth to convert to screen space radius
	s.add_uniform_vector3d("center", center);
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

void add_depth_of_field(float focus_depth, float dof_val) {

	bind_depth_buffer(1); // tu_id=1
	shader_t s;

	for (unsigned dim = 0; dim < 2; ++dim) {
		bind_frame_buffer_RGB();
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("depth_utils.part+depth_of_field");
		s.begin_shader();
		setup_depth_tex(s, 1);
		s.add_uniform_int  ("dim_val",     (dim ? 1 : 0));
		s.add_uniform_float("focus_depth", focus_depth);
		s.add_uniform_float("dof_val",     dof_val);
		fill_screen_white_and_end_shader(s);
		color_buffer_frame = 0; // reset to invalidate buffer and force recreation of texture for second pass
	}
}

void add_2d_blur() { // faster than add_color_only_effect("screen_space_blur")

	shader_t s;

	for (unsigned dim = 0; dim < 2; ++dim) {
		bind_frame_buffer_RGB();
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("screen_space_blur_1d");
		s.begin_shader();
		s.add_uniform_int("frame_buffer_tex", 0);
		set_xy_step(s);
		s.add_uniform_int("dim_val", (dim ? 1 : 0));
		fill_screen_white_and_end_shader(s);
		color_buffer_frame = 0; // reset to invalidate buffer and force recreation of texture for second pass
	}
}

void add_bloom() {

	bind_frame_buffer_RGB();
	shader_t s;

	for (unsigned dim = 0; dim < 2; ++dim) {
		if (dim) {bind_frame_buffer_RGB();}
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("postproc_bloom");
		s.begin_shader();
		s.add_uniform_int("frame_buffer_tex", 0);
		set_xy_step(s);
		s.add_uniform_int("dim_val", (dim ? 1 : 0));
		fill_screen_white_and_end_shader(s);
		color_buffer_frame = 0; // reset to invalidate buffer and force recreation of texture for second pass
	}
}

void add_2d_bloom() {

	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("postproc_bloom_2d");
	s.begin_shader();
	s.add_uniform_int("frame_buffer_tex", 0);
	set_xy_step(s);
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer and force recreation of texture for second pass
}

void add_postproc_underwater_fog(float atten_scale, float max_uw_dist, float mud_amt, float algae_amt) {

	bind_depth_buffer(1); // tu_id=1
	shader_t s;
	bind_frame_buffer_RGB();
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("depth_utils.part+postproc_water_fog");
	s.begin_shader();
	s.add_uniform_float("max_uw_dist", max_uw_dist);
	setup_depth_tex(s, 1);
	setup_shader_underwater_atten(s, atten_scale, mud_amt, algae_amt);
	fill_screen_white_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}

void apply_player_underwater_effect(colorRGBA const &color_mod=WHITE, float intensity=1.0) {
	add_2d_blur();
	if (player_is_drowning())                 {add_color_only_effect("drunken_wave", 1.00*intensity, 1.0, 0.0, color_mod);}
	else if (world_mode == WMODE_INF_TERRAIN) {add_color_only_effect("drunken_wave", 0.12*intensity, 1.6, 1.6, color_mod);} // reduced but faster effect
}

void run_postproc_effects() {

	bool const enable_ssao = 0;
	point const camera(get_camera_pos());
	bool const camera_underwater(world_mode != WMODE_UNIVERSE && is_underwater(camera));
	float const drunkenness(get_player_drunkenness());
	int index(-1);
	static xform_matrix prev_mvm, prev_pjm; // previous frame's matrices, for use with motion blur, etc.
	static bool prev_mat_valid(0);
	if (enable_ssao && (display_mode & 0x20)) {add_ssao();} // key '6'
	
	if (cur_explosion_sphere.radius > 0.0 && camera_pdu.sphere_visible_test(cur_explosion_sphere.pos, cur_explosion_sphere.radius)) {
		if (dist_less_than(camera, cur_explosion_sphere.pos, max(40.0*cur_explosion_sphere.radius, 8.0*CAMERA_RADIUS))) { // close/large on the screen
			if (coll_pt_vis_test(camera, cur_explosion_sphere.pos, 0.0, index, camera_coll_id, 1, 1)) {
				add_sphere_refract_effect(cur_explosion_sphere, cur_explosion_weight);
			}
		}
	}
	if (drunkenness > 0.5) { // at least slightly drunk
		if (drunkenness > 1.5) {add_2d_blur();} // very drunk
		if (drunkenness > 1.0) {add_color_only_effect("double_vision", 0.5f*(drunkenness - 1.0f));} // moderately drunk
		add_color_only_effect("drunken_wave", 1.0f*(min(drunkenness, 1.25f) - 0.5f));
	}
	else if (camera_underwater) {
		apply_player_underwater_effect();
	}
	else {
		float const dist_to_fire(sqrt(dist_to_fire_sq)), fire_max_dist(4.0*CAMERA_RADIUS);
		float const dist_to_lava((water_is_lava && world_mode == WMODE_INF_TERRAIN) ? (camera.z - water_plane_z) : 0.0f), lava_max_dist(16.0*CAMERA_RADIUS);
		if      (dist_to_fire > 0.0 && dist_to_fire < fire_max_dist) {add_color_only_effect("heat_waves", (fire_max_dist - dist_to_fire)/fire_max_dist);}
		else if (dist_to_lava > 0.0 && dist_to_lava < lava_max_dist) {add_color_only_effect("heat_waves", (lava_max_dist - dist_to_lava)/lava_max_dist);}
	}
	if (display_mode & 0x80) { // DOF, key '8'
		point const pos2(camera + cview_dir*FAR_CLIP);
		point cpos(pos2);
		vector3d cnorm; // unused
		int cindex(-1); // unused
		float focus_depth(FAR_CLIP);
		if (check_coll_line_exact(camera, pos2, cpos, cnorm, cindex, 0.0, camera_coll_id)) {focus_depth = p2p_dist(camera, cpos);}
		if (line_intersect_mesh(camera, cpos, cpos)) {focus_depth = p2p_dist(camera, cpos);}
		float const dof_val(0.04*FAR_CLIP);
		add_depth_of_field(focus_depth, dof_val);
	}
	if (show_fog && world_mode == WMODE_GROUND && !camera_underwater && !is_rain_enabled()) {add_god_rays();}
	
	if (!enable_ssao && (display_mode & 0x20) && !camera_underwater) { // add bloom last, key '6'
		if (world_mode != WMODE_INF_TERRAIN) {add_bloom();}
		else if (have_buildings() && is_night()) {add_2d_bloom();} // allow bloom for building windows at night in TT mode
	}
	if (enable_postproc_recolor) {add_color_only_effect("recolor", 0.0);} // add recolor at the very end
	if (vignette_color.A > 0.0 ) {add_vignette(vignette_color);}

	if (0 && !prev_mat_valid) { // capture matrices from this frame for use with next frame (if needed in the future)
		prev_mvm = fgGetMVM();
		prev_pjm = fgGetPJM();
		prev_mat_valid = 1;
	}
}

