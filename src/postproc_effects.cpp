// 3D World - Full screen postprocessing effects
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "function_registry.h"
#include "gl_ext_arb.h"
#include "shaders.h"


extern unsigned depth_tid, frame_buffer_RGB_tid;
extern int frame_counter, display_mode, show_fog, camera_coll_id, window_width, window_height;
extern float NEAR_CLIP, FAR_CLIP;
extern colorRGBA sun_color;


void bind_depth_buffer(bool force_regen=0) {

	static int prev_frame_counter(-1);
	if (force_regen || frame_counter != prev_frame_counter) {depth_buffer_to_texture(depth_tid);} // depth texture is not valid for this frame
	prev_frame_counter = frame_counter;
	assert(depth_tid >= 0);
	bind_2d_texture(depth_tid);
}

void bind_frame_buffer_RGB(bool force_regen=0) {
	
	static int prev_frame_counter(-1);
	if (force_regen || frame_counter != prev_frame_counter) {frame_buffer_RGB_to_texture(frame_buffer_RGB_tid);} // FB RGB texture is not valid for this frame
	prev_frame_counter = frame_counter;
	assert(frame_buffer_RGB_tid >= 0);
	bind_2d_texture(frame_buffer_RGB_tid);
}

void draw_ortho_screen_space_quad() {

	// setup matrices
	fgMatrixMode(FG_PROJECTION);
	fgPushIdentityMatrix();
	fgOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	fgMatrixMode(FG_MODELVIEW);
	fgPushIdentityMatrix();

	enable_blend();
	glDisable(GL_DEPTH_TEST);
	draw_tquad(1.0, 1.0, 0.0);
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

void draw_white_quad_and_end_shader(shader_t &s) {

	s.set_cur_color(WHITE);
	draw_ortho_screen_space_quad();
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
	draw_white_quad_and_end_shader(s);
}

void add_ssao() {

	bind_depth_buffer();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("depth_utils.part+screen_space_ao");
	s.begin_shader();
	setup_depth_tex(s, 0);
	draw_white_quad_and_end_shader(s);
}

void add_color_blur() {

	bind_frame_buffer_RGB();
	shader_t s;
	s.set_vert_shader("no_lighting_tex_coord");
	s.set_frag_shader("screen_space_blur");
	s.begin_shader();
	s.add_uniform_int("frame_buffer_tex", 0);
	set_xy_step(s);
	draw_white_quad_and_end_shader(s);
}

void add_depth_of_field(float focus_depth, float dof_val) {

	set_active_texture(1);
	bind_depth_buffer();
	set_active_texture(0);
	bind_frame_buffer_RGB();
	shader_t s;

	for (unsigned dim = 0; dim < 2; ++dim) {
		if (dim) {bind_frame_buffer_RGB(1);} // force recreation of texture for second pass
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("depth_utils.part+depth_of_field");
		s.begin_shader();
		setup_depth_tex(s, 1);
		s.add_uniform_float("dim_val",     (dim ? 1.0 : 0.0));
		s.add_uniform_float("focus_depth", focus_depth);
		s.add_uniform_float("dof_val",     dof_val);
		draw_white_quad_and_end_shader(s);
	}
}

void run_postproc_effects() {

	point const camera(get_camera_pos());
	if (show_fog && world_mode == WMODE_GROUND) {add_god_rays();}
	if (display_mode & 0x20) {add_ssao();}
	if (world_mode != WMODE_UNIVERSE && is_underwater(camera)) {add_color_blur();}
	
	if (display_mode & 0x80) {
		point const pos2(camera + cview_dir*FAR_CLIP);
		point cpos(pos2);
		vector3d cnorm; // unused
		int cindex(-1), xpos, ypos; // unused
		float focus_depth(FAR_CLIP), zval;
		if (check_coll_line_exact(camera, pos2, cpos, cnorm, cindex, 0.0, camera_coll_id)) {focus_depth = p2p_dist(camera, cpos);}
		if (line_intersect_mesh(camera, cpos, xpos, ypos, zval)) {focus_depth = p2p_dist(camera, (camera + (cpos - camera)*(zval - camera.z)/(cpos.z - camera.z)));}
		float const dof_val(0.04*FAR_CLIP);
		add_depth_of_field(focus_depth, dof_val);
	}
}

