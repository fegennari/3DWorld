// 3D World - Drawing Code
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "dynamic_particle.h"
#include "physics_objects.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "draw_utils.h"
#include "transform_obj.h"
#include <glm/vec4.hpp>


bool const DYNAMIC_SMOKE_SHADOWS = 1; // slower, but looks nice
bool const MIN_PARTICLE_FILL     = 1;
unsigned const MAX_CFILTERS      = 10;
float const NDIV_SCALE           = 1.6;
float const CLOUD_WIND_SPEED     = 0.00015;


struct sky_pos_orient {

	point center;
	float radius, radius_inv, dx, dy;
	sky_pos_orient(point const &c, float r, float dx_, float dy_)
		: center(c), radius(r), radius_inv(1.0/radius), dx(dx_), dy(dy_) {assert(radius > 0.0);}
};


bool enable_clip_plane_z(0);
unsigned depth_tid(0), frame_buffer_RGB_tid(0);
float sun_radius(0.0), moon_radius(0.0), earth_radius(0.0), brightness(1.0), clip_plane_z(0.0);
colorRGB cur_ambient(BLACK), cur_diffuse(BLACK);
point sun_pos, moon_pos;
gl_light_params_t gl_light_params[MAX_SHADER_LIGHTS];
point const earth_pos(-22.0, -12.0, 31.0);
sky_pos_orient cur_spo(point(0,0,0),1,0,0);
vector3d up_norm(plus_z);
vector<camera_filter> cfilters;
pt_line_drawer bubble_pld;

extern bool have_sun, using_lightmap, has_dl_sources, has_spotlights, has_line_lights, smoke_exists, two_sided_lighting, tree_indir_lighting;
extern bool group_back_face_cull, have_indir_smoke_tex, combined_gu, enable_depth_clamp, dynamic_smap_bias, volume_lighting, dl_smap_enabled, underwater;
extern bool enable_gamma_correct, smoke_dlights;
extern int is_cloudy, iticks, frame_counter, display_mode, show_fog, use_smoke_for_fog, num_groups, xoff, yoff;
extern int window_width, window_height, game_mode, draw_model, camera_mode, DISABLE_WATER, animate2, camera_coll_id;
extern unsigned smoke_tid, dl_tid, create_voxel_landscape, enabled_lights, reflection_tid, scene_smap_vbo_invalid, sky_zval_tid;
extern float zmin, light_factor, fticks, perspective_fovy, perspective_nclip, cobj_z_bias;
extern float temperature, atmosphere, zbottom, indir_vert_offset, rain_wetness, snow_cov_amt, NEAR_CLIP, FAR_CLIP;
extern point light_pos, mesh_origin, flow_source, surface_pos;
extern vector3d wind;
extern colorRGB const_indir_color, ambient_lighting_scale;
extern colorRGBA bkg_color, sun_color, sunlight_color, base_cloud_color, cur_fog_color;
extern lightning l_strike;
extern vector<spark_t> sparks;
extern vector<star> stars;
extern vector<beam3d> beams;
extern set<unsigned> moving_cobjs;
extern obj_group obj_groups[];
extern coll_obj_group coll_objects;
extern cobj_draw_groups cdraw_groups;
extern obj_type object_types[];
extern obj_vector_t<bubble> bubbles;
extern obj_vector_t<particle_cloud> part_clouds;
extern cloud_manager_t cloud_manager;
extern obj_vector_t<fire> fires;
extern obj_vector_t<decal_obj> decals;
extern water_particle_manager water_part_man;
extern physics_particle_manager explosion_part_man;
extern cube_t cur_smoke_bb;
extern vector<portal> portals;
extern vector<obj_draw_group> obj_draw_groups;
extern reflect_plane_selector reflect_planes;

void create_dlight_volumes();
void create_sky_vis_zval_texture(unsigned &tid);
cube_t get_all_models_bcube(bool only_reflective=0);


void set_fill_mode() {
	glPolygonMode(GL_FRONT_AND_BACK, ((draw_model == 0) ? GL_FILL : GL_LINE));
}

void ensure_filled_polygons() {
	if (draw_model != 0) {glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);} // always filled
}

void reset_fill_mode() {
	if (draw_model != 0) {set_fill_mode();}
}

int get_universe_ambient_light(bool for_universe_draw) {return (for_universe_draw ? 1 : 3);}


void gl_light_params_t::set_pos(point const &p, float w) {
	
	// Note: it might seem like we can skip the update of eye_space_pos when pos is unchanged;
	// however, eye_space_pos depends on the MVM, which can change independently (and does for ship engines, etc.)
	//if (p == pos && w == pos_w) return 0; // already set to this value
	pos   = p;
	pos_w = w;
	glm::vec4 const v(fgGetMVM() * glm::vec4(p.x, p.y, p.z, w));
	eye_space_pos.assign(v.x, v.y, v.z); // v.w ignored?
}


void set_light_ds_color(int light, colorRGBA const &diffuse, shader_t *shader) {

	assert(light >= 0 && light < MAX_SHADER_LIGHTS);
	gl_light_params[light].set_ds(diffuse);
	if (shader) {shader->upload_light_source(light, 0x0C);}
}

void set_light_a_color(int light, colorRGBA const &ambient, shader_t *shader) {

	assert(light >= 0 && light < MAX_SHADER_LIGHTS);
	gl_light_params[light].set_a(ambient);
	if (shader) {shader->upload_light_source(light, 0x02);}
}

void set_light_colors(int light, colorRGBA const &ambient, colorRGBA const &diffuse, shader_t *shader) {

	set_light_ds_color(light, diffuse, shader);
	set_light_a_color (light, ambient, shader);
}

void set_colors_and_enable_light(int light, colorRGBA const &ambient, colorRGBA const &diffuse, shader_t *shader) {

	enable_light(light);
	set_light_colors(light, ambient, diffuse, shader);
}

void clear_colors_and_disable_light(int light, shader_t *shader) {

	enable_light(light); // enable temporarily so that we can update the shader colors for shaders that don't check the enabled state
	assert(light >= 0 && light < MAX_SHADER_LIGHTS);
	set_light_colors(light, BLACK, BLACK, shader);
	disable_light(light);
}

void set_gl_light_pos(int light, point const &pos, float w, shader_t *shader) {

	assert(light >= 0 && light < MAX_SHADER_LIGHTS);
	gl_light_params[light].set_pos(pos, w);
	if (shader) {shader->upload_light_source(light, 0x01);}
}

void setup_gl_light_atten(int light, float c_a, float l_a, float q_a, shader_t *shader) {

	gl_light_params[light].set_atten(c_a, l_a, q_a);
	if (shader) {shader->upload_light_source(light, 0x10);}
}


// metal material properties - low diffuse color but very high specular of a similar color to diffuse
// Note: diffuse color should be zero, but is set to nonzero (half albedo) due to lack of env maps, which contribute multi-directional specular reflections
void set_silver_material(shader_t &shader, float alpha) {
	shader.set_cur_color(colorRGBA(WHITE*0.5, alpha));
	shader.set_specular_color(colorRGBA(0.8, 0.8, 0.8)*2.0, 50.0);
}
void set_gold_material(shader_t &shader, float alpha) {
	shader.set_cur_color(colorRGBA(GOLD*0.5, alpha));
	shader.set_specular_color(colorRGBA(0.9, 0.6, 0.1)*2.0, 40.0);
}
void set_copper_material(shader_t &shader, float alpha) {
	shader.set_cur_color(colorRGBA(COPPER_C*0.5, alpha));
	shader.set_specular_color(colorRGBA(0.8, 0.4, 0.25)*2.0, 30.0);
}
void set_brass_material(shader_t &shader, float alpha) {
	shader.set_cur_color(colorRGBA(BRASS_C*0.5, alpha));
	shader.set_specular_color(colorRGBA(0.85, 0.85, 0.25)*2.0, 25.0);
}


void draw_camera_weapon(bool want_has_trans) {

	if (!game_mode || weap_has_transparent(CAMERA_ID) != want_has_trans) return;
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1);
	draw_weapon_in_hand(-1, s);
	s.end_shader();
}


bool is_light_enabled(int l) {assert(l < MAX_SHADER_LIGHTS); return ((enabled_lights & (1<<l)) != 0);}
void enable_light    (int l) {assert(l < MAX_SHADER_LIGHTS); enabled_lights |=  (1<<l);}
void disable_light   (int l) {assert(l < MAX_SHADER_LIGHTS); enabled_lights &= ~(1<<l);}


void calc_cur_ambient_diffuse() {

	unsigned ncomp(0);
	cur_ambient = cur_diffuse = BLACK;
	if (l_strike.enabled) {cur_ambient += LITN_C*0.25;}

	for (unsigned i = 0; i < 2; ++i) { // sun, moon
		if (!is_light_enabled(i)) continue;
		cur_ambient += gl_light_params[i].ambient;
		cur_diffuse += gl_light_params[i].diffuse;
		++ncomp;
	}
	if (ncomp > 0) {
		float const cscale(0.5 + 0.5/ncomp);
		cur_ambient  = cur_ambient.modulate_with(ambient_lighting_scale) * cscale;
		cur_diffuse *= cscale;
	}
}


bool set_dlights_booleans(shader_t &s, bool enable, int shader_type, bool no_dl_smap) {

	if (!enable) {s.set_prefix("#define NO_DYNAMIC_LIGHTS", shader_type);} // if we're not even enabling dlights
	bool const dl_en(enable && dl_tid > 0 && has_dl_sources);
	
	if (dl_en) {
		if (has_spotlights)  {s.set_prefix("#define HAS_SPOTLIGHTS",  shader_type);}
		if (has_line_lights) {s.set_prefix("#define HAS_LINE_LIGHTS", shader_type);}
		if (dl_smap_enabled && !no_dl_smap) {s.set_prefix("#define HAS_DLIGHT_SMAP", shader_type);}
	}
	s.set_prefix(make_shader_bool_prefix("enable_dlights", dl_en), shader_type);
	return dl_en;
}


void common_shader_block_pre(shader_t &s, bool &dlights, bool &use_shadow_map, bool &indir_lighting, float min_alpha, bool no_dl_smap, bool use_wet_mask=0) {

	bool const hemi_lighting(!have_indir_smoke_tex);
	use_shadow_map &= shadow_map_enabled();
	indir_lighting &= have_indir_smoke_tex;
	dlights        &= (dl_tid > 0 && has_dl_sources);
	s.check_for_fog_disabled();
	if (enable_gamma_correct) {s.set_prefix("#define ENABLE_GAMMA_CORRECTION", 1);} // FS
	if (min_alpha == 0.0    ) {s.set_prefix("#define NO_ALPHA_TEST", 1);} // FS
	if (use_shadow_map && dynamic_smap_bias) {s.set_prefix("#define DYNAMIC_SMAP_BIAS", 1);} // FS
	s.set_prefix(make_shader_bool_prefix("indir_lighting", indir_lighting), 1); // FS
	s.set_prefix(make_shader_bool_prefix("hemi_lighting",  hemi_lighting),  1); // FS
	s.set_prefix(make_shader_bool_prefix("use_shadow_map", use_shadow_map), 1); // FS
	s.set_prefix(make_shader_bool_prefix("use_water_coverage", use_wet_mask), 1); // FS
	set_dlights_booleans(s, dlights, 1, no_dl_smap); // FS
}


void set_indir_color(shader_t &s) {
	colorRGB const indir_color((have_indir_smoke_tex && world_mode == WMODE_GROUND) ? colorRGB(0.0, 0.0, 0.0) : const_indir_color);
	s.add_uniform_color("const_indir_color", indir_color);
}

void set_indir_lighting_block(shader_t &s, bool use_smoke, bool use_indir) {

	s.setup_scene_bounds();
	if ((use_smoke || use_indir) && smoke_tid) {set_3d_texture_as_current(smoke_tid, 1);}
	s.add_uniform_int("smoke_and_indir_tex", 1);
	s.add_uniform_float("half_dxy", HALF_DXY);
	s.add_uniform_float("indir_vert_offset", indir_vert_offset);
	s.add_uniform_float("ambient_scale", (use_indir ? 0.0 : 1.0)); // ambient handled by indirect lighting in the shader
	set_indir_color(s);

	// hemispherical lighting
	s.add_uniform_color("sky_color", colorRGB(bkg_color));
	select_multitex(LANDSCAPE_TEX, 12, 1);
	s.add_uniform_int("ground_tex", 12);
}


void common_shader_block_post(shader_t &s, bool dlights, bool use_shadow_map, bool use_smoke, bool use_indir, float min_alpha, bool enable_dlights_smap=1) {

	s.setup_fog_scale(); // fog scale for the case where smoke is disabled
	if (dlights) {setup_dlight_textures(s, enable_dlights_smap);}
	set_indir_lighting_block(s, use_smoke, use_indir);
	s.add_uniform_int("tex0", 0);
	s.add_uniform_float("min_alpha", min_alpha);
	// the z plane bias is somewhat of a hack, set experimentally; maybe should be one pixel in world space?
	if (enable_clip_plane_z) {s.add_uniform_float("clip_plane_z", (clip_plane_z - 0.005));} // 10.0*cobj_z_bias? 0.2/window_height?
	if (use_shadow_map && shadow_map_enabled() && world_mode == WMODE_GROUND) {set_smap_shader_for_all_lights(s, cobj_z_bias);}
	set_active_texture(0);
	s.clear_specular();
}


float const SMOKE_NOISE_MAG = 0.8;
bool use_smoke_noise() {return (use_smoke_for_fog == 2 && SMOKE_NOISE_MAG);}
bool is_smoke_in_use() {return (smoke_exists || use_smoke_for_fog);}


void set_smoke_shader_prefixes(shader_t &s, int use_texgen, bool keep_alpha, bool direct_lighting,
	bool smoke_enabled, int has_lt_atten, bool use_smap, int use_bmap, bool use_spec_map, bool use_mvm, bool use_tsl)
{
	s.set_int_prefix("use_texgen", use_texgen, 0); // VS
	s.set_prefix(make_shader_bool_prefix("keep_alpha",          keep_alpha),      1); // FS
	s.set_prefix(make_shader_bool_prefix("direct_lighting",     direct_lighting), 1); // FS
	s.set_prefix(make_shader_bool_prefix("do_cube_lt_atten",    ((has_lt_atten & 1) != 0)), 1); // FS
	s.set_prefix(make_shader_bool_prefix("do_sphere_lt_atten",  ((has_lt_atten & 2) != 0)), 1); // FS
	s.set_prefix(make_shader_bool_prefix("two_sided_lighting",  use_tsl), 1); // FS
	s.set_prefix(make_shader_bool_prefix("use_fg_ViewMatrix",   use_mvm), 0); // VS
	s.set_prefix(make_shader_bool_prefix("enable_clip_plane_z", enable_clip_plane_z), 1); // FS
	if (use_spec_map) {s.set_prefix("#define USE_SPEC_MAP",   1);} // FS
	if (use_smap && shadow_map_enabled()) {s.set_prefix("#define USE_SHADOW_MAP", 1);} // FS
	
	if (smoke_enabled) {
		// Note: dynamic_smoke_shadows applies to light0 only
		// Note: dynamic_smoke_shadows still uses the visible smoke bbox, so if you can't see smoke it won't cast a shadow
		for (unsigned d = 0; d < 2; ++d) { // VS/FS
			if (DYNAMIC_SMOKE_SHADOWS) {s.set_prefix("#define DYNAMIC_SMOKE_SHADOWS", d);}
			s.set_prefix("#define SMOKE_ENABLED", d);
		}
		if (volume_lighting && use_smap) {s.set_prefix("#define SMOKE_SHADOW_MAP", 1);} // FS
		if (smoke_dlights) {s.set_prefix("#define SMOKE_DLIGHTS", 1);} // FS - TESTING

		if (use_smoke_noise()) {
			s.set_prefix("#define SMOKE_NOISE",   1); // FS
			s.set_prefix("#define NUM_OCTAVES 3", 1); // FS
			//s.set_prefix("#define RIDGED_NOISE",  1); // FS
		}
	}
	if (use_bmap) {
		for (unsigned i = 0; i < 2; ++i) {
			s.set_prefix("#define USE_BUMP_MAP", i); // VS/FS
			if (use_bmap == 2) {s.set_prefix("#define USE_TANGENT_VECTOR", i);} // VS/FS
		}
		s.set_prefix("#define USE_BUMP_MAP_INDIR", 1); // FS
		s.set_prefix("#define USE_BUMP_MAP_DL",    1); // FS
	}
	s.setup_enabled_lights(3, 2); // FS; sun, moon, and lightning
}


float setup_underwater_fog(shader_t &s, int shader_type) {
	float water_depth(0.0);
	bool const underwater((world_mode == WMODE_GROUND) ? is_underwater(get_camera_pos(), 1, &water_depth) : 0);
	s.set_prefix(make_shader_bool_prefix("underwater", underwater), shader_type);
	return water_depth;
}

unsigned get_sky_zval_texture() {

	if (sky_zval_tid == 0) {create_sky_vis_zval_texture(sky_zval_tid);}
	assert(sky_zval_tid != 0);
	return sky_zval_tid;
}

void invalidate_snow_coverage() {free_texture(sky_zval_tid);}


// texture units used: 0: object texture, 1: smoke/indir lighting texture, 2-4 dynamic lighting, 5: bump map, 6-7: shadow map,
//                     8: specular map, 9: depth map, 10: burn mask/sky_zval, 11: noise, 12: ground texture, 13: depth, 14: reflection, 15: ripples, 16-31: dlight shadow maps
// use_texgen: 0 = use texture coords, 1 = use standard texture gen matrix, 2 = use custom shader tex0_s/tex0_t, 3 = use vertex id for texture
// use_bmap  : 0 = none, 1 = auto generate tangent vector, 2 = tangent vector in vertex attribute
// is_outside: 0 = inside, 1 = outside, 2 = use snow coverage mask
void setup_smoke_shaders(shader_t &s, float min_alpha, int use_texgen, bool keep_alpha, bool indir_lighting, bool direct_lighting, bool dlights, bool smoke_en,
	int has_lt_atten, bool use_smap, int use_bmap, bool use_spec_map, bool use_mvm, bool force_tsl, float burn_tex_scale, float triplanar_texture_scale,
	bool use_depth_trans, bool enable_reflections, int is_outside, bool enable_rain_snow)
{
	bool const triplanar_tex(triplanar_texture_scale != 0.0);
	bool const use_burn_mask(burn_tex_scale > 0.0);
	bool const is_wet(is_ground_wet() && !use_burn_mask), is_snowy(enable_rain_snow && is_ground_snowy() && !use_burn_mask), use_wet_mask(is_wet && is_outside == 2);
	bool const enable_puddles(enable_rain_snow && is_wet && !is_rain_enabled()); // enable puddles when the ground is wet but it's not raining
	smoke_en &= (have_indir_smoke_tex && smoke_tid > 0 && is_smoke_in_use());
	if (use_burn_mask     ) {s.set_prefix("#define APPLY_BURN_MASK",        1);} // FS
	if (triplanar_tex     ) {s.set_prefix("#define TRIPLANAR_TEXTURE",      1);} // FS
	if (use_depth_trans   ) {s.set_prefix("#define USE_DEPTH_TRANSPARENCY", 1);} // FS
	if (enable_reflections) {s.set_prefix("#define ENABLE_REFLECTIONS",     1);} // FS
	if (enable_puddles    ) {s.set_prefix("#define ENABLE_PUDDLES",         1);} // FS
	if (is_snowy          ) {s.set_prefix("#define ENABLE_SNOW_COVERAGE",   1);} // FS
	float const water_depth(setup_underwater_fog(s, 1)); // FS
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, min_alpha, 0, use_wet_mask);
	set_smoke_shader_prefixes(s, use_texgen, keep_alpha, direct_lighting, smoke_en, has_lt_atten, use_smap, use_bmap, use_spec_map, use_mvm, force_tsl);
	s.set_vert_shader("texture_gen.part+bump_map.part+leaf_wind.part+no_lt_texgen_smoke");
	string fstr("linear_fog.part+bump_map.part+spec_map.part+ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+line_clip.part*+indir_lighting.part+black_body_burn.part+");
	if (smoke_en && use_smoke_noise()) {fstr += "perlin_clouds_3d.part*+";}
	if (enable_reflections) {fstr += "water_ripples.part+";}
	if (triplanar_tex) {fstr += "triplanar_texture.part+";}
	if (use_depth_trans) {fstr += "depth_utils.part+";}
	s.set_frag_shader(fstr + "textured_with_smoke");
	s.begin_shader();
	s.add_uniform_float("water_depth", water_depth);

	if (use_texgen == 2) {
		s.register_attrib_name("tex0_s", TEX0_S_ATTR);
		s.register_attrib_name("tex0_t", TEX0_T_ATTR);
	}
	if (use_bmap == 2) {s.register_attrib_name("tangent", TANGENT_ATTR);}
	if (use_bmap     ) {s.add_uniform_int("bump_map", 5);}
	if (use_spec_map ) {s.add_uniform_int("spec_map", 8);}
	if (triplanar_tex) {s.add_uniform_float("tex_scale", triplanar_texture_scale);}
	common_shader_block_post(s, dlights, use_smap, smoke_en, indir_lighting, min_alpha);
	float step_delta_scale((use_smoke_for_fog || get_smoke_at_pos(get_camera_pos())) ? 1.0 : 2.0);
	s.add_uniform_float("step_delta_shadow", step_delta_scale*HALF_DXY);
	if (volume_lighting && is_light_enabled(0)) {step_delta_scale *= 0.2f;} // 5 steps per texel for sun light on smoke volume
	s.add_uniform_float("step_delta", step_delta_scale*HALF_DXY);
	if (use_mvm) {upload_mvm_to_shader(s, "fg_ViewMatrix");}
	
	if (smoke_en) {
		cube_t const smoke_bb(get_scene_bounds());
		s.add_uniform_float_array("smoke_bb", (use_smoke_for_fog ? &smoke_bb.d[0][0] : &cur_smoke_bb.d[0][0]), 6);
		if (DYNAMIC_SMOKE_SHADOWS) {s.add_uniform_vector3d("sun_pos", get_sun_pos());}
		s.add_uniform_color("smoke_color",     (use_smoke_for_fog ? colorRGB(cur_fog_color) : colorRGB(GRAY)));
		s.add_uniform_float("smoke_const_add", ((use_smoke_for_fog == 1) ? 0.25 : 0.0));

		if (use_smoke_noise()) {
			set_3d_texture_as_current(get_noise_tex_3d(64, 1), 11); // grayscale noise
			s.add_uniform_int("cloud_noise_tex", 11);
			s.add_uniform_float("smoke_noise_mag", SMOKE_NOISE_MAG);
			s.add_uniform_float("noise_scale", 0.45);
			static vector3d fog_time(zero_vector);
			static int update_frame(0);
			if (animate2 && frame_counter > update_frame) {fog_time += 0.001*fticks*wind; update_frame = frame_counter;} // fog moves with the wind (once per frame)
			s.add_uniform_vector3d("fog_time", fog_time);
		}
	}
	if (use_burn_mask) {
		s.add_uniform_float("burn_tex_scale", burn_tex_scale);
		s.add_uniform_float("burn_offset", -1.0); // starts disabled
		s.add_uniform_int("burn_mask", 10);
		select_multitex(DISINT_TEX, 10); // PLASMA_TEX?
	}
	if (enable_reflections) {
		s.add_uniform_int("reflection_tex", 14);
		select_multitex(RIPPLE_MAP_TEX, 15);
		s.add_uniform_int("ripple_tex", 15);
		static float ripple_time(0.0);
		static int update_frame(0);
		if (animate2 && frame_counter > update_frame) {ripple_time += fticks; update_frame = frame_counter;} // once per frame
		cube_t const &bcube(reflect_planes.get_selected());
		s.add_uniform_float("ripple_time",        ripple_time);
		s.add_uniform_float("rain_intensity",     get_rain_intensity());
		s.add_uniform_float("reflect_plane_zbot", bcube.d[2][0]);
		s.add_uniform_float("reflect_plane_ztop", bcube.d[2][1]);
	}
	if (enable_puddles) {
		set_3d_texture_as_current(get_noise_tex_3d(64, 1), 11); // grayscale noise
		s.add_uniform_int("wet_noise_tex", 11);
	}
	if (use_wet_mask || is_snowy) {
		bind_texture_tu(get_sky_zval_texture(), 10);
		s.add_uniform_int("sky_zval_tex", 10);
	}
	// need to handle wet/outside vs. dry/inside surfaces differently, so the caller must either set is_outside properly or override wet and snow values
	s.add_uniform_float("wet_effect",   (is_outside ? rain_wetness : 0.0)); // only enable when drawing cobjs?
	s.add_uniform_float("reflectivity", (enable_reflections ?  1.0 : 0.0));
	s.add_uniform_float("snow_cov_amt", snow_cov_amt); // Note: no longer depends on is_outside
}


void set_tree_branch_shader(shader_t &s, bool direct_lighting, bool dlights, bool use_smap) {

	float const water_depth(setup_underwater_fog(s, 1)); // FS
	bool indir_lighting(direct_lighting && tree_indir_lighting);
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, 0.0, 1); // no_dl_smap=1
	set_smoke_shader_prefixes(s, 0, 0, direct_lighting, 0, 0, use_smap, 0, 0, 0, 0);
	s.set_vert_shader("texture_gen.part+bump_map.part+leaf_wind.part+no_lt_texgen_smoke");
	s.set_frag_shader("linear_fog.part+bump_map.part+ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+line_clip.part*+indir_lighting.part+textured_with_smoke");
	s.begin_shader();
	s.add_uniform_float("water_depth", water_depth);
	common_shader_block_post(s, dlights, use_smap, 0, indir_lighting, 0.0, 0); // no dlights smap
	check_gl_error(400);
}


// texture units used: 0,8,15: object texture, 1: indir lighting texture, 2-4: dynamic lighting, 5: 3D noise texture, 6-7: shadow map, 9: AO texture, 10: voxel shadow texture, 11: normal map texture, 12: ground texture
void setup_procedural_shaders(shader_t &s, float min_alpha, bool indir_lighting, bool dlights, bool use_smap, bool use_bmap,
	bool use_noise_tex, bool z_top_test, float tex_scale, float noise_scale, float tex_mix_saturate)
{
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, min_alpha, 0);
	
	if (use_bmap) {
		// FIXME: only looks correct with sun/moon lighting - dynamic and indirect lighting doesn't work with triplanar bump mapping for some reason, but it's too slow anyway
		s.set_prefix("#define USE_BUMP_MAP",    1); // FS
		s.set_prefix("#define BUMP_MAP_CUSTOM", 1); // FS
		//s.set_prefix("#define USE_BUMP_MAP_INDIR", 1); // FS
		//s.set_prefix("#define USE_BUMP_MAP_DL",    1); // FS
	}
	s.set_prefix(make_shader_bool_prefix("use_noise_tex", use_noise_tex), 1); // FS
	s.set_prefix(make_shader_bool_prefix("z_top_test",    z_top_test),    1); // FS
	s.setup_enabled_lights(2, 2); // FS; only 2, but could be up to 8 later
	s.set_vert_shader("procedural_gen");
	s.set_frag_shader("linear_fog.part+bump_map.part+ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+triplanar_texture.part+procedural_texture.part+indir_lighting.part+voxel_texture.part+triplanar_bump_map.part+procedural_gen");
	s.begin_shader();
	common_shader_block_post(s, dlights, use_smap, 0, indir_lighting, min_alpha);
	s.add_uniform_int("tex1",    8);
	s.add_uniform_int("tex_top", 15); // not used in all cases
	s.add_uniform_float("tex_scale", tex_scale);

	if (use_noise_tex) {
		s.add_uniform_int("noise_tex", 5); // does this need an enable option?
		s.add_uniform_float("noise_scale", noise_scale);
		s.add_uniform_float("tex_mix_saturate", tex_mix_saturate);
	}
	if (use_bmap) {
		select_multitex(ROCK_NORMAL_TEX, 11, 1);
		s.add_uniform_int("bump_map", 11);
		s.add_uniform_float("bump_tex_scale", 4.0);
	}
}


void setup_object_render_data() {

	RESET_TIME;
	bool const TIMETEST(0);
	calc_cur_ambient_diffuse();
	create_dlight_volumes();
	distribute_smoke();
	if (TIMETEST) {PRINT_TIME("1 Distribute Smoke");}
	upload_smoke_indir_texture();
	if (TIMETEST) {PRINT_TIME("2 Upload Smoke");}
	add_dynamic_lights_ground(); // and create dlights shadow maps
	if (TIMETEST) {PRINT_TIME("3 Add Dlights");}
	cube_t const dlight_bounds(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, get_zval_min(), get_zval_max());
	upload_dlights_textures(dlight_bounds); // get_scene_bounds()
	if (TIMETEST) {PRINT_TIME("4 Dlights Textures");}
	get_occluders();
	if (TIMETEST) {PRINT_TIME("5 Get Occluders");}
	//scene_smap_vbo_invalid = 0; // needs to be after dlights update
}


bool  enable_reflection_plane() {return ((display_mode & 0x10) && !reflect_planes.empty());}
bool  use_reflection_plane   () {return (enable_reflection_plane() && reflect_planes.enabled() && get_camera_pos().z > reflect_planes.get_selected().d[2][0]);}
float get_reflection_plane   () {return reflect_planes.get_refl_plane();}

bool use_reflect_plane_for_cobj(coll_obj const &c) {
	if (c.type != COLL_CUBE && (c.type != COLL_CYLINDER || (c.cp.surfs & 1))) return 0;
	if (!c.is_wet() && c.cp.spec_color.get_luminance() < 0.25) return 0;
	cube_t const &bc(reflect_planes.get_selected());
	return (c.intersects(bc) && c.d[2][1] <= bc.d[2][1] && camera_pdu.cube_visible(c));
}

void proc_refl_bcube(cube_t const &c, cube_t &bcube, float &min_camera_dist, bool &bcube_set) {

	point const camera(get_camera_pos());
	float const dist(p2p_dist(camera, c.closest_pt(camera)));
		
	if (bcube_set) {
		bcube.union_with_cube(c);
		min_camera_dist = min(dist, min_camera_dist);
	}
	else {
		bcube = c;
		bcube_set = 1;
		min_camera_dist = dist;
	}
}

bool get_reflection_plane_bounds(cube_t &bcube, float &min_camera_dist) {

	reflect_planes.select_best_reflection_plane();
	if (!use_reflection_plane()) return 0;
	bool bcube_set(0);
	point const camera(get_camera_pos());

	for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
		unsigned cix(*i);
		coll_obj const &c(coll_objects.get_cobj(cix));
		if (c.no_draw() || c.d[2][1] >= camera.z || !use_reflect_plane_for_cobj(c)) continue;
		cube_t cc(c);
		cc.d[2][0] = cc.d[2][1]; // shrink to top surface only
		proc_refl_bcube(cc, bcube, min_camera_dist, bcube_set);
	}
	cube_t const models_refl_bcube(get_all_models_bcube(1));
	
	if (!models_refl_bcube.is_zero_area() && models_refl_bcube.intersects(reflect_planes.get_selected())) {
		proc_refl_bcube(models_refl_bcube, bcube, min_camera_dist, bcube_set);
	}
	return bcube_set;
}


void end_group(int &last_group_id) {

	if (last_group_id < 0) return;
	assert((unsigned)last_group_id < obj_draw_groups.size());
	obj_draw_groups[last_group_id].end_render();
	if (group_back_face_cull) glDisable(GL_CULL_FACE);
	last_group_id = -1;
}

coll_obj const &get_draw_cobj(unsigned index) {
	if (index >= coll_objects.size()) {return cdraw_groups.get_cobj(index - coll_objects.size());}
	return coll_objects.get_cobj(index);
}

void setup_cobj_shader(shader_t &s, bool has_lt_atten, bool enable_normal_maps, int use_texgen, bool enable_reflections) {
	// Note: pass in 3 when has_lt_atten to enable sphere atten
	setup_smoke_shaders(s, 0.0, use_texgen, 0, 1, 1, 1, 1, has_lt_atten, 1, enable_normal_maps, 0, (use_texgen == 0), two_sided_lighting, 0.0, 0.0, 0, enable_reflections);
}

void draw_cobjs_group(vector<unsigned> const &cobjs, cobj_draw_buffer &cdb, shader_t &s, int use_texgen, bool use_normal_map, bool use_reflect_tex) {

	if (cobjs.empty()) return;
	setup_cobj_shader(s, 0, use_normal_map, use_texgen, use_reflect_tex); // no lt_atten
	if (use_reflect_tex) {bind_texture_tu(reflection_tid, 14);}
	cdb.clear();
	// we use generated tangent and binormal vectors, with the binormal scale set to either 1.0 or -1.0 depending on texture coordinate system and y-inverting
	float bump_b_scale(0.0);
	int nm_tid(-2), last_tid(-2), last_group_id(-1); // Note: use -2 as unset tid so that it differs from "no texture" of -1

	// Note: could stable_sort normal_map_cobjs by normal_map tid, but the normal map is already part of the layer sorting,
	// so the normal maps are probably already grouped together
	for (auto i = cobjs.begin(); i != cobjs.end(); ++i) {
		coll_obj const &c(get_draw_cobj(*i));

		if (use_normal_map) {
			float const bbs(c.cp.negate_nm_bns() ? -1.0 : 1.0);

			if (c.cp.normal_map != nm_tid) { // normal map change
				cdb.flush();
				nm_tid = c.cp.normal_map;
				select_multitex(nm_tid, 5);
			}
			if (bbs != bump_b_scale) {
				cdb.flush();
				bump_b_scale = bbs;
				s.add_uniform_float("bump_b_scale", bump_b_scale);
			}
		}
		unsigned cix(*i);
		c.draw_cobj(cix, last_tid, last_group_id, s, cdb);
		assert(cix == *i); // should not have been modified
		//if (use_reflect_tex) {bind_2d_texture(reflection_tid);} // overwrite texture binding (FIXME: temporary)
	}
	cdb.flush();
	s.clear_specular(); // may be unnecessary
	if (use_normal_map) {s.add_uniform_float("bump_b_scale", -1.0);} // reset, may be unnecessary
	s.end_shader();
}

typedef vector<pair<float, int> > vect_sorted_ix;

bool check_big_occluder(coll_obj const &c, unsigned cix, vect_sorted_ix &out) { // Note: increases CPU time but decreases GPU time

	if (!c.is_big_occluder() || c.group_id >= 0) return 0;
	float const dist_sq(distance_to_camera_sq(c.get_center_pt()));
	if (c.get_area() < 0.005*dist_sq) return 0;
	if (!camera_pdu.cube_visible(c)) return 1; // skip this cobj - return 1, but don't add as an occluder
	out.push_back(make_pair(dist_sq, cix));
	return 1;
}

bool draw_or_add_cobj(unsigned cix, bool reflection_pass, bool use_ref_plane, vect_sorted_ix large_cobjs[2], vect_sorted_ix &draw_last,
	vector<unsigned> &normal_map_cobjs, vector<unsigned> &tex_coord_cobjs, vector<unsigned> &tex_coord_nm_cobjs,
	vector<unsigned> &reflect_cobjs, vector<unsigned> &reflect_cobjs_nm)
{
	coll_obj const &c(get_draw_cobj(cix));
	bool const cube_poly((c.cp.flags & COBJ_WAS_CUBE) != 0);
	bool const use_tex_coords((cube_poly && c.cp.tid >= 0) || ((c.is_cylinder() || c.type == COLL_SPHERE || c.type == COLL_CAPSULE) && c.cp.tscale == 0.0));

	// Note: only texgen cube/vert cylinder top surfaces support reflections
	if (!use_tex_coords && use_ref_plane && use_reflect_plane_for_cobj(c)) {
		assert(c.group_id < 0);
		if (reflection_pass) return 0; // the reflection surface is not drawn in the reflection pass (receiver only)
		((c.cp.normal_map >= 0) ? reflect_cobjs_nm : reflect_cobjs).push_back(cix);
		return 0;
	}
	if (c.cp.normal_map >= 0) { // common case
		assert(c.group_id < 0);
		assert(!c.is_semi_trans());
		if (!use_tex_coords && check_big_occluder(c, cix, large_cobjs[1])) return 0;
		(use_tex_coords ? tex_coord_nm_cobjs : normal_map_cobjs).push_back(cix);
		return 0;
	}
	if (use_tex_coords) { // uncommon case (typically movable objects)
		assert(c.group_id < 0);
		//assert(!c.is_semi_trans()); // too strong?
		if (camera_pdu.cube_visible(c)) {tex_coord_cobjs.push_back(cix);}
		return 0;
	}
	if (check_big_occluder(c, cix, large_cobjs[0])) return 0;

	if (c.is_semi_trans()) { // slow when polygons are grouped
		float dist(distance_to_camera(c.get_center_pt()));
		if (c.type == COLL_SPHERE) {dist -= c.radius;} // distance to surface closest to the camera
		else if (c.type == COLL_CYLINDER || c.type == COLL_CYLINDER_ROT || c.type == COLL_CAPSULE) { // approx distance to surface closest to the camera
			dist -= min(0.5*(c.radius + c.radius2), 0.5*p2p_dist(c.points[0], c.points[1]));
		}
		draw_last.push_back(make_pair(-dist, cix)); // negative distance
		return 0;
	}
	return 1;
}

// should always have draw_solid enabled on the first call for each frame
void draw_coll_surfaces(bool draw_trans, bool reflection_pass) {

	//RESET_TIME;
	static vect_sorted_ix draw_last;
	if (coll_objects.empty() || coll_objects.drawn_ids.empty() || world_mode != WMODE_GROUND) return;
	if (draw_trans && draw_last.empty() && (!is_smoke_in_use() || portals.empty())) return; // nothing transparent to draw
	// Note: in draw_solid mode, we could call get_shadow_triangle_verts() on occluders to do a depth pre-pass here, but that doesn't seem to be more efficient
	bool const has_lt_atten(draw_trans && coll_objects.has_lt_atten);
	bool const use_ref_plane(reflection_pass || (reflection_tid > 0 && use_reflection_plane()));
	float const ref_plane_z(use_ref_plane ? get_reflection_plane() : 0.0);
	shader_t s;
	setup_cobj_shader(s, has_lt_atten, 0, 2, 0);
	int last_tid(-2), last_group_id(-1);
	static cobj_draw_buffer cdb;
	cdb.clear();
	vector<unsigned> normal_map_cobjs, tex_coord_cobjs, tex_coord_nm_cobjs, reflect_cobjs, reflect_cobjs_nm;
	//if (enable_clip_plane_z) {glEnable(GL_CLIP_DISTANCE0);}
	
	if (!draw_trans) { // draw solid
		vect_sorted_ix large_cobjs[2]; // {normal_map=0, normal_map=1}
		draw_last.clear();

		for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
			unsigned cix(*i);
			coll_obj const &c(coll_objects.get_cobj(cix));
			assert(c.cp.draw);
			if (c.no_draw()) continue; // can still get here sometimes
			if (reflection_pass && c.d[2][1] < ref_plane_z) continue; // below the reflection plane (approximate) (optimization)
			assert(c.id == (int)cix); // should always be equal

			if (c.dgroup_id >= 0) {
				if (!c.is_cobj_visible()) continue; // VFC/occlusion culling
				vector<unsigned> const &group_cids(cdraw_groups.get_draw_group(c.dgroup_id, c));
				//assert(!group_cids.empty()); // too strong?

				for (auto j = group_cids.begin(); j != group_cids.end(); ++j) {
					unsigned const ix(*j + coll_objects.size()); // map to a range that doesn't overlap coll_objects

					if (draw_or_add_cobj(ix, reflection_pass, use_ref_plane, large_cobjs, draw_last, normal_map_cobjs, tex_coord_cobjs, tex_coord_nm_cobjs, reflect_cobjs, reflect_cobjs_nm)) {
						unsigned ix2(ix);
						cdraw_groups.get_cobj(*j).draw_cobj(ix2, last_tid, last_group_id, s, cdb); // Note: ix should not be modified
						assert(ix2 == ix); // should not have changed
					}
				}
				continue; // don't draw c itself (only if group is nonempty?)
			}
			if (draw_or_add_cobj(cix, reflection_pass, use_ref_plane, large_cobjs, draw_last, normal_map_cobjs, tex_coord_cobjs, tex_coord_nm_cobjs, reflect_cobjs, reflect_cobjs_nm)) {
				c.draw_cobj(cix, last_tid, last_group_id, s, cdb); // i may not be valid after this call
				
				if (cix != *i) {
					assert(cix > *i);
					i = std::lower_bound(i, coll_objects.drawn_ids.end(), cix);
				}
			}
		} // for i
		//cout << "cobjs: " << coll_objects.drawn_ids.size() << ", nm: " << normal_map_cobjs.size() << ", tc: " << tex_coord_cobjs.size() << ", tcnm: " << tex_coord_nm_cobjs.size() << ", lg0: " << large_cobjs[0].size() << ", lg1: " << large_cobjs[1].size() << ", trans: " << draw_last.size() << endl;
		end_group(last_group_id);
		cdb.flush();
		for (unsigned d = 0; d < 2; ++d) {sort(large_cobjs[d].begin(), large_cobjs[d].end());} // sort front to back for early z culling

		for (auto i = large_cobjs[0].begin(); i != large_cobjs[0].end(); ++i) { // not normal mapped
			unsigned cix(i->second);
			get_draw_cobj(cix).draw_cobj(cix, last_tid, last_group_id, s, cdb);
		}
		cdb.flush();
		for (auto i = large_cobjs[1].begin(); i != large_cobjs[1].end(); ++i) {normal_map_cobjs.push_back(i->second);} // normal mapped
	} // end draw solid
	else { // draw transparent
		if (is_smoke_in_use()) {
			for (unsigned i = 0; i < portals.size(); ++i) {
				if (!portals[i].is_visible()) continue;
				draw_last.push_back(make_pair(-distance_to_camera(portals[i].get_center_pt()), -(int)(i+1)));
			}
		}
		sort(draw_last.begin(), draw_last.end()); // sort back to front for alpha blending
		enable_blend();
		lt_atten_manager_t lt_atten_manager(s);
		if (has_lt_atten) {lt_atten_manager.enable();}
		vector<vert_wrap_t> portal_verts;
		bool in_portal(0);
		// disable depth writing so that other partially transparent objects drawn in the wrong alpha order (fires, smoke, decals, particles)
		// aren't discarded by the depth buffer, even though they won't be properly alpha blended, assuming alpha value is small (windows, portals, etc.)
		glDepthMask(GL_FALSE);

		for (unsigned i = 0; i < draw_last.size(); ++i) {
			int const ix(draw_last[i].second);

			if (ix < 0) { // portal
				end_group(last_group_id);
				cdb.flush();
				lt_atten_manager.next_object(0.0, 1.0);
				if (!in_portal) {portal::pre_draw(portal_verts); in_portal = 1;}
				unsigned const pix(-(ix+1));
				assert(pix < portals.size());
				portals[pix].draw(portal_verts);
			}
			else { // cobj
				if (in_portal) {portal::post_draw(portal_verts); in_portal = 0;}
				unsigned cix(ix);
				coll_obj const &c(get_draw_cobj(cix));
				cdb.on_new_obj_layer(c.cp);
				bool using_lt_atten(0);
				
				if (has_lt_atten) { // we only support cubes and spheres for now (Note: may not be compatible with groups)
					if (c.type == COLL_CUBE) {
						lt_atten_manager.next_cube(c.cp.light_atten, c.cp.refract_ix, c);
						using_lt_atten = (c.cp.light_atten > 0.0);
					}
					else if (c.type == COLL_SPHERE) {
						lt_atten_manager.next_sphere(c.cp.light_atten, c.cp.refract_ix, c.points[0], c.radius);
						using_lt_atten = (c.cp.light_atten > 0.0);
					}
					else {
						lt_atten_manager.next_object(0.0, c.cp.refract_ix); // reset
					}
				}
				c.draw_cobj(cix, last_tid, last_group_id, s, cdb);
				if (using_lt_atten) {cdb.flush();} // must flush because ulocs[2] is per-cube
				assert((int)cix == ix); // should not have changed
			}
		} // for i
		if (in_portal) {portal::post_draw(portal_verts);}
		end_group(last_group_id);
		cdb.flush();
		disable_blend();
		draw_last.resize(0);
		glDepthMask(GL_TRUE); // re-enable depth writing
	} // end draw_trans
	s.clear_specular(); // may be unnecessary
	s.end_shader();
	draw_cobjs_group(normal_map_cobjs,   cdb, s, 2, 1, 0);
	draw_cobjs_group(tex_coord_nm_cobjs, cdb, s, 0, 1, 0);
	draw_cobjs_group(tex_coord_cobjs,    cdb, s, 0, 0, 0);
	draw_cobjs_group(reflect_cobjs,      cdb, s, 2, 0, 1);
	draw_cobjs_group(reflect_cobjs_nm,   cdb, s, 2, 1, 1);
	//if (enable_clip_plane_z) {glDisable(GL_CLIP_DISTANCE0);}
	//if (draw_solid) PRINT_TIME("Final Draw");
}


bool portal::is_visible() const {

	cube_t bcube(pts, 4);
	if (normal != zero_vector && dot_product_ptv(normal, get_camera_pos(), bcube.get_cube_center()) < 0.0) return 0; // back facing
	if (!camera_pdu.cube_visible(bcube)) return 0;
	if ((display_mode & 0x08) && cobj_contained(get_camera_pos(), pts, 4, -1)) return 0;
	return 1;
}

void portal::pre_draw(vector<vert_wrap_t> &verts) {

	select_texture(WHITE_TEX);
	ALPHA0.set_for_cur_shader();
	assert(verts.empty());
}

void portal::post_draw(vector<vert_wrap_t> &verts) {

	draw_verts(verts, GL_TRIANGLES);
	verts.clear();
}

void portal::draw(vector<vert_wrap_t> &verts) const {

	for (unsigned i = 0; i < 6; ++i) {verts.push_back(pts[quad_to_tris_ixs[i]]);}
};


void draw_stars(float alpha) {

	if (alpha <= 0.0) return;
	colorRGBA color(BLACK), bkg;
	UNROLL_3X(bkg[i_] = (1.0 - alpha)*bkg_color[i_];)
	point const xlate((camera_mode == 1) ? surface_pos : all_zeros);
	enable_blend();
	glDisable(GL_DEPTH_TEST);
	point_sprite_drawer psd;
	psd.reserve_pts(stars.size());
	rand_gen_t rgen;

	for (unsigned i = 0; i < stars.size(); ++i) {
		if ((rgen.rand()%400) == 0) continue; // flicker out

		for (unsigned j = 0; j < 3; ++j) {
			float const c(stars[i].color[j]*stars[i].intensity);
			color[j] = ((alpha >= 1.0) ? c : (alpha*c + bkg[j]));
		}
		psd.add_pt(vert_color((stars[i].pos + xlate), color));
	}
	psd.draw(BLUR_TEX, 2.0); // draw with points of size 2 pixels
	glEnable(GL_DEPTH_TEST);
	disable_blend();
}


void draw_sun() {

	point const pos(get_sun_pos());
	if (!have_sun || !sphere_in_camera_view(pos, sun_radius, 2)) return;
	colorRGBA color(attenuate_sun_color(SUN_C));
	apply_red_sky(color);
	draw_single_colored_sphere(pos, sun_radius, N_SPHERE_DIV, color);
}


void draw_moon() {

	if (world_mode == WMODE_GROUND && show_fog) return; // don't draw when there is fog
	point const pos(get_moon_pos());
	if (!sphere_in_camera_view(pos, moon_radius, 2)) return;
	colorRGBA const ambient(0.05, 0.05, 0.05, 1.0), diffuse(1.0*have_sun, 1.0*have_sun, 1.0*have_sun, 1.0);
	set_colors_and_enable_light(4, ambient, diffuse);
	set_gl_light_pos(4, get_sun_pos(), 0.0);
	shader_t s;
	s.set_vert_shader("moon_draw");
	s.set_frag_shader("ads_lighting.part*+moon_draw");
	s.begin_shader();
	s.add_uniform_int("tex0", 0);
	s.set_cur_color(attenuate_sun_color(WHITE));
	select_texture(MOON_TEX);
	draw_cube_mapped_sphere(pos, moon_radius, N_SPHERE_DIV/2, 1);
	s.end_shader();
	disable_light(4);
	float const star_alpha(get_star_alpha());

	if (star_alpha < 1.0) { // fade moon into background color when the sun comes up
		enable_blend();
		draw_single_colored_sphere(pos, 1.1*moon_radius, N_SPHERE_DIV, colorRGBA(bkg_color, (1.0 - star_alpha)));
		disable_blend();
	}
}


// for some reason the texture is backwards, so we mirrored the image of the earth
void draw_earth() {

	if (show_fog) return; // don't draw when there is fog
	point pos(mesh_origin + earth_pos);
	if (camera_mode == 1) {pos += surface_pos;}
	static float rot_angle(0.0);

	if (sphere_in_camera_view(pos, earth_radius, 2)) {
		colorRGBA const color(attenuate_sun_color(WHITE));
		shader_t s;
		s.begin_simple_textured_shader(0.0, 1, 0, &color);
		select_texture(EARTH_TEX);
		fgPushMatrix();
		translate_to(pos);
		fgRotate(67.0, 0.6, 0.8, 0.0);
		fgRotate(rot_angle, 0.0, 0.0, 1.0);
		fgRotate(180.0, 1.0, 0.0, 0.0);
		draw_sphere_vbo(all_zeros, earth_radius, N_SPHERE_DIV, 1); // Note: texture is pre-distorted for sphere mapping
		fgPopMatrix();
		s.end_shader();
	}
	if (animate2) {rot_angle += 0.2*fticks;}
}


void apply_red_sky(colorRGBA &color) {

	if (light_factor > 0.45 && light_factor < 0.55) { // red sky at night/morning
		float const redness(1.0 - 20.0*fabs(light_factor - 0.5));
		color.R = min(1.0f, (1.0f + 0.8f*redness)*color.R);
		color.G = max(0.0f, (1.0f - 0.2f*redness)*color.G);
		color.B = max(0.0f, (1.0f - 0.5f*redness)*color.B);
	}
}


colorRGBA get_cloud_color() {

	colorRGBA color(brightness, brightness, min(1.0, 1.2*brightness), atmosphere); // more blue when cloudy/rainy
	apply_red_sky(color);
	return color.modulate_with(base_cloud_color);
}


void get_avg_sky_color(colorRGBA &avg_color) {
	blend_color(avg_color, colorRGBA(get_cloud_color(), 1.0), bkg_color, 0.5, 1);
}


void set_cloud_intersection_shader(shader_t &s) {

	s.add_uniform_vector3d("sun_pos",  sun_pos);
	s.add_uniform_vector3d("moon_pos", moon_pos);
	s.add_uniform_vector3d("sphere_center", cur_spo.center);
	s.add_uniform_float   ("sphere_radius", cur_spo.radius);
	s.add_uniform_float   ("atmosphere", atmosphere);
	s.add_uniform_float   ("dx", cur_spo.dx);
	s.add_uniform_float   ("dy", cur_spo.dy);
	s.add_uniform_int("cloud_tex", 8);
	select_multitex(CLOUD_TEX, 8, 1);
}


float get_cloud_density(point const &pt, vector3d const &dir) { // optimize?

	if (atmosphere == 0.0) return 0.0;
	point lsint;
	if (!line_sphere_int(-dir, pt, cur_spo.center, cur_spo.radius, lsint, 0)) return 0.0; // shouldn't get here?
	vector3d const vdir(lsint - cur_spo.center);
	return atmosphere*get_texture_component(CLOUD_TEX, (vdir.x*cur_spo.radius_inv + cur_spo.dx), (vdir.y*cur_spo.radius_inv + cur_spo.dy), 3); // cloud alpha
}


void draw_sky(bool camera_side, bool no_update) {

	if (atmosphere < 0.01) return; // no atmosphere
	float radius(0.55*(FAR_CLIP+X_SCENE_SIZE));
	point center((camera_mode == 1) ? surface_pos : mesh_origin);
	center.z -= 0.727*radius;
	if ((distance_to_camera(center) > radius) != camera_side) return;
	colorRGBA const cloud_color(get_cloud_color());
	static float sky_rot_xy[2] = {0.0, 0.0}; // x, y
	float const wmag(sqrt(wind.x*wind.x + wind.y*wind.y));

	if (!no_update && wmag > TOLERANCE) {
		for (unsigned d = 0; d < 2; ++d) {sky_rot_xy[d] -= fticks*CLOUD_WIND_SPEED*(wmag + 0.5*WIND_ADJUST)*wind[d]/wmag;}
	}
	cur_spo = sky_pos_orient(center, radius, sky_rot_xy[0], sky_rot_xy[1]);
	int const light_id(4);
	enable_blend();

	if (have_sun && light_factor > 0.4) { // sun lighting of clouds
		colorRGBA horizon_color;
		float const blend_val(atmosphere*CLIP_TO_01(2.0f*(1.0f - get_star_alpha())));
		blend_color(horizon_color, WHITE, ALPHA0, blend_val, 1);
		horizon_color.alpha *= 0.5;
		apply_red_sky(horizon_color);
		horizon_color = horizon_color.modulate_with(base_cloud_color);
		shader_t s;
		s.begin_simple_textured_shader(0.0, 0, 0, &horizon_color);
		select_texture(GRADIENT_TEX);
		draw_sphere_vbo(center, 1.05*radius, N_SPHERE_DIV, 1); // draw horizon
		s.end_shader();

		point lpos(get_sun_pos()), lsint;
		vector3d const sun_v((get_camera_pos() - lpos).get_norm());
		if (line_sphere_int(sun_v, lpos, center, radius, lsint, 1)) {lpos = lsint;}
		colorRGBA const ambient(sun_color*0.5);
		set_colors_and_enable_light(light_id, ambient, sun_color);
		set_gl_light_pos(light_id, lpos, 1.0); // w=1.0 - point light source
		setup_gl_light_atten(light_id, 0.0, 0.01, 0.01);
	}
	if (enable_depth_clamp) {glDisable(GL_DEPTH_CLAMP);}
	glDepthMask(GL_FALSE); // disable depth writing
	shader_t s;
	s.setup_enabled_lights(5, 1); // sun, moon, and L4 VS lighting (L2 and L3 are set but unused)
	s.set_vert_shader("ads_lighting.part*+texture_gen.part+cloud_sphere");
	s.set_frag_shader("simple_texture");
	s.begin_shader();
	s.add_uniform_float("min_alpha", 0.0);
	s.add_uniform_int("tex0", 0);
	s.clear_specular();
	s.set_cur_color(cloud_color);
	select_texture(CLOUD_TEX);
	// change S and T parameters to map sky texture into the x/y plane with translation based on wind/rot
	setup_texgen(1.0, 1.0, sky_rot_xy[0], sky_rot_xy[1], 0.0, s, 0); // Note: in post-transformed sphere space
	draw_sphere_vbo(center, radius, (3*N_SPHERE_DIV)/2, 1);
	s.end_shader();
	disable_blend();
	disable_light(light_id);
	glDepthMask(GL_TRUE); // enable depth writing
	if (enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);}
}


void compute_brightness() {

	brightness = 0.8 + 0.2*light_factor;
	if (!have_sun) {brightness *= 0.25;}
	if (is_cloudy) {brightness *= 0.5;}
	float const sun_bright(0.5 + 0.5*max(0.0f, sun_pos.z/sun_pos.mag()));
	float const moon_bright(combined_gu ? 0.1 : 0.3*(0.5 + 0.5*max(0.0f, moon_pos.z/moon_pos.mag())));
	if      (light_factor >= 0.6) {brightness *= sun_bright;}
	else if (light_factor <= 0.4) {brightness *= moon_bright;}
	else {brightness *= 5.0*((light_factor - 0.4)*sun_bright + (0.6 - light_factor)*moon_bright);}
	brightness = min(0.99f, max(0.0f, brightness));
}


template<typename S, typename T> void get_draw_order(vector<T> const &objs, vector<S> &order) {

	point const camera(get_camera_pos());
	
	for (unsigned i = 0; i < objs.size(); ++i) {
		if (!objs[i].status) continue;
		point const pos(objs[i].get_pos());
		if (sphere_in_camera_view(pos, objs[i].radius, 0)) {order.push_back(make_pair(-p2p_dist_sq(pos, camera), i));}
	}
	sort(order.begin(), order.end()); // sort back to front
}


void bubble::draw(bool set_liquid_color) const {

	assert(status);
	colorRGBA color2(color);
	if (set_liquid_color) {select_liquid_color(color2, pos);}
	float const point_dia(NDIV_SCALE*window_width*radius/distance_to_camera(pos));

	if (point_dia < 4.0) {
		bubble_pld.add_pt(pos, (get_camera_pos() - pos), color2);
	}
	else {
		color2.set_for_cur_shader();
		int const ndiv(max(4, min(16, int(4.0*sqrt(point_dia)))));
		draw_sphere_vbo(pos, radius, ndiv, 0);
	}
}


order_vect_t particle_cloud::order;


void particle_cloud::draw(quad_batch_draw &qbd) const {

	assert(status);
	colorRGBA color(base_color);
	color.A *= density;
	
	if (is_fire()) {
		color.G *= get_rscale();
	}
	else {
		color *= (no_lighting ? 1.0 : brightness)*(0.5*(1.0 - darkness));
	}
	if (parts.empty()) {
		if (status && sphere_in_camera_view(pos, radius, 0)) {draw_part(pos, radius, color, qbd);}
	}
	else {
		order.resize(0);
		render_parts.resize(parts.size());

		for (unsigned i = 0; i < parts.size(); ++i) {
			render_parts[i].pos    = pos + parts[i].pos*radius;
			render_parts[i].radius = parts[i].radius*radius;
			render_parts[i].status = parts[i].status;
		}
		get_draw_order(render_parts, order);
		
		for (unsigned j = 0; j < order.size(); ++j) {
			unsigned const i(order[j].second);
			assert(i < render_parts.size());
			draw_part(render_parts[i].pos, render_parts[i].radius, color, qbd);
		}
	}
}


void particle_cloud::draw_part(point const &p, float r, colorRGBA c, quad_batch_draw &qbd) const {

	point const camera(get_camera_pos());
	if (dist_less_than(camera, p, max(NEAR_CLIP, 4.0f*r))) return; // too close to the camera

	if (!no_lighting && !is_fire()) { // fire has its own emissive lighting
		point const lpos(get_light_pos());
		bool const outside_scene(p.z > czmax || !is_over_mesh(p));
		bool known_coll(0);
		static int last_cid(-1);
		if (!outside_scene && last_cid >= 0) {known_coll = coll_objects.get_cobj(last_cid).line_intersect(p, lpos);}

		if (outside_scene || (!known_coll && !check_coll_line(p, lpos, last_cid, -1, 1, 1))) { // not shadowed (slow, especially for lots of smoke near trees)
			// Note: This can be moved into a shader, but the performance and quality improvement might not be significant
			vector3d const dir((p - get_camera_pos()).get_norm());
			float const dp(dot_product_ptv(dir, p, lpos));
			float rad, dist, t;
			colorRGBA const cloud_color(get_cloud_color());
			blend_color(c, cloud_color, c, 0.15, 0); // 15% ambient lighting (transmitted/scattered)
			if (dp > 0.0) {blend_color(c, cloud_color, c, 0.1*dp/p2p_dist(p, lpos), 0);} // 10% diffuse lighting (directional)

			if (dp < 0.0 && have_sun && line_intersect_sphere(p, dir, sun_pos, 6*sun_radius, rad, dist, t)) {
				float const mult(1.0 - max(0.0f, (rad - sun_radius)/(5*sun_radius)));
				blend_color(c, SUN_C, c, 0.75*mult, 0); // 75% direct sun lighting
			}
		}
		get_indir_light(c, p); // could move outside of the parts loop if too slow
	}
	if (red_only) c.G = c.B = 0.0; // for special luminosity cloud texture rendering
	// Note: Can disable smoke volume integration for close smoke, but very close smoke (< 1 grid unit) is infrequent
	qbd.add_billboard(p, camera, up_vector, c, 4.0*r, 4.0*r, tex_range_t(), (MIN_PARTICLE_FILL && !no_lighting)); // use quads for clouds
}


void fire::draw(quad_batch_draw &qbd, int &last_in_smoke) const {

	assert(status);
	point const pos2(pos + point(0.0, 0.0, 2.0*radius));
	int const in_smoke((get_smoke_at_pos(get_camera_pos()) || get_smoke_at_pos(pos2)) != 0.0);

	if (in_smoke != last_in_smoke) {
		qbd.draw_and_clear();
		if (in_smoke) {set_std_blend_mode();} else {set_additive_blend_mode();}
		last_in_smoke = in_smoke;
	}
	qbd.add_animated_billboard(pos2, get_camera_pos(), up_vector, WHITE, 4.0*radius, 4.0*radius, (time&15)/16.0);
}


bool decal_obj::draw(quad_batch_draw &qbd) const {

	assert(status);
	point const cur_pos(get_pos());
	bool const back_facing(dot_product_ptv(orient, cur_pos, get_camera_pos()) > 0.0);
	if (back_facing && (cid < 0 || !coll_objects[cid].cp.is_glass())) return 0; // back face culling
	float const alpha_val(get_alpha());
	if (!dist_less_than(cur_pos, get_camera_pos(), max(window_width, window_height)*radius*alpha_val)) return 0; // distance culling
	vector3d upv(orient.y, orient.z, orient.x); // swap the xyz values to get an orthogonal vector
	rotate_vector3d(orient, rot_angle, upv);
	// move slightly away from the object to blend properly with cracks
	point const viewer(cur_pos + (back_facing ? -orient : orient));
	qbd.add_billboard((cur_pos + DECAL_OFFSET*orient), viewer, upv, colorRGBA(color, alpha_val), radius, radius, tex_range);
	return 1;
}

// drawn a line of blood running down a vertical surface
void decal_obj::maybe_draw_blood_trail(line_tquad_draw_t &blood_tqd) const {

	assert(status);
	if (orient.z != 0.0)  return; // not a vertical surface
	if (color != BLOOD_C) return; // not blood
	if (cid < 0)          return; // no cobj (can this happen?)
	int const hashval(int(1000.0*pos.z));
	if (hashval & 1)      return; // only half of them have blood
	coll_obj const &cobj(coll_objects.get_cobj(cid));
	if (cobj.type != COLL_CUBE) return; // only cubes for now
	point const cur_pos(get_pos()), camera(get_camera_pos());
	if (!cobj.cp.is_glass() && dot_product_ptv(orient, cur_pos, camera) > 0.0) return; // back face culling
	float const alpha_val(get_alpha()), view_thresh(max(window_width, window_height)*radius*alpha_val);
	if (!dist_less_than(cur_pos, get_camera_pos(), 0.7*view_thresh)) return; // distance culling
	float const dz_max(min(min((cur_pos.z - cobj.d[2][0]), 100.0f*radius), 2.0f*CAMERA_RADIUS)), dz(((hashval&63)+1)/64.0 * dz_max);
	if (dz < radius) return; // too small
	if (!sphere_in_camera_view(cur_pos-vector3d(0,0, 0.5*dz), radius+0.5*dz, 0)) return; // VFC
	colorRGBA const c(color, alpha_val);
	float const width(0.15*radius);
	point const start(cur_pos - vector3d(0,0, 1.0*radius)), end(cur_pos - vector3d(0,0, dz)), end2(end - vector3d(0,0, 0.15*radius));
	blood_tqd.add_line_as_tris(start, end, width, width, c, c); // trail
	if (dist_less_than(cur_pos, camera, 0.25*view_thresh)) {blood_tqd.add_line_as_tris(cur_pos, start, 0.6*width, width, colorRGBA(c, 0.0), c);} // transition
	if (dist_less_than(end,     camera, 0.10*view_thresh)) {blood_tqd.add_line_as_tris(end, end2, width, 0.0, c, c);} // end - circle?
}


template<typename T, typename ARG> void draw_objects(vector<T> const &objs, ARG &arg) {

	order_vect_t order;
	get_draw_order(objs, order);

	for (unsigned i = 0; i < order.size(); ++i) {
		assert(order[i].second < objs.size());
		objs[order[i].second].draw(arg);
	}
}


void draw_bubbles() {

	if (!bubbles.any_active()) return;
	shader_t s;
	s.begin_untextured_lit_glcolor_shader();
	glEnable(GL_CULL_FACE);
	enable_blend();
	bool const set_liquid_color(world_mode == WMODE_GROUND);
	begin_sphere_draw(0);
	draw_objects(bubbles, set_liquid_color);
	end_sphere_draw();
	bubble_pld.draw_and_clear();
	disable_blend();
	glDisable(GL_CULL_FACE);
	s.end_shader();
}


void draw_part_clouds(vector<particle_cloud> const &pc, int tid) {

	enable_flares(tid);
	//select_multitex(CLOUD_TEX, 1);
	static quad_batch_draw qbd;
	draw_objects(pc, qbd);
	qbd.draw_and_clear(); // color will be set per object
	disable_flares();
	//set_active_texture(0);
}


void physics_particle_manager::draw(float radius, int tid) const {

	if (parts.empty()) return;
	point const camera(get_camera_pos());
	enable_blend();
	point_sprite_drawer_norm_sized psd;
	psd.reserve_pts(parts.size());

	for (unsigned i = 0; i < parts.size(); ++i) {
		//psd.add_pt(vert_norm_color(parts[i].p, (camera - parts[i].p).get_norm(), parts[i].c.c), radius); // normal faces camera
		psd.add_pt(sized_vert_t<vert_norm_color>(vert_norm_color(parts[i].p, (camera - parts[i].p).get_norm(), parts[i].c.c), radius)); // normal faces camera
	}
	if (tid >= 0) {psd.sort_back_to_front();} // if we have an alpha texture, sort back to front
	psd.draw(tid, 0.0, 1); // draw with lighting
	disable_blend();
}

void water_particle_manager::draw() const {
	physics_particle_manager::draw(0.5*object_types[DROPLET].radius, BLUR_CENT_TEX); // constant value, half the size of regular droplets
}


struct crack_point {

	point pos, orig_pos;
	int cid, face, time;
	float alpha;
	colorRGBA color;
	
	crack_point() {}
	crack_point(point const &pos_, point const &opos, int cid_, int face_, int time_, float alpha_, colorRGBA const &color_)
		: pos(pos_), orig_pos(opos), cid(cid_), face(face_), time(time_), alpha(alpha_), color(color_) {}
	
	bool operator<(crack_point const &c) const {
		if (cid  != c.cid ) return (cid  < c.cid );
		if (face != c.face) return (face < c.face);
		return (c.time < time); // max time first
	}
};


struct ray2d {

	point2d<float> pts[2];

	ray2d() {}
	ray2d(float x1, float y1, float x2, float y2) {pts[0].x = x1; pts[0].y = y1; pts[1].x = x2; pts[1].y = y2;}
};


void create_and_draw_cracks(quad_batch_draw &qbd) { // adds to beams

	vector<crack_point> cpts;  // static?
	vector<ray2d> crack_lines; // static?
	int last_cobj(-1);
	bool skip_cobj(0);
	point const camera(get_camera_pos());

	for (vector<decal_obj>::const_iterator i = decals.begin(); i != decals.end(); ++i) {
		if (i->status == 0 || !i->is_glass || i->cid < 0)    continue;
		if (i->cid == last_cobj && skip_cobj)                continue;
		if (moving_cobjs.find(i->cid) != moving_cobjs.end()) continue;
		point const pos(i->get_pos());
		if (!dist_less_than(camera, pos, 2000*i->radius))    continue; // too far away
		coll_obj const &cobj(coll_objects.get_cobj(i->cid));
		skip_cobj = (cobj.status != COLL_STATIC || cobj.type != COLL_CUBE || !camera_pdu.cube_visible(cobj) || cobj.is_occluded_from_camera());
		last_cobj = i->cid;
		if (skip_cobj) continue;
		int const face(cobj.closest_face(pos));
		if ((face >> 1) != get_min_dim(cobj)) continue; // only draw cracks in the plane of the min dimension/thickness
		cpts.push_back(crack_point(pos, i->pos, i->cid, face, i->time, i->get_alpha(), i->color));
	}
	stable_sort(cpts.begin(), cpts.end());

	for (unsigned i = 0; i < cpts.size();) {
		unsigned const s(i);
		for (++i; i < cpts.size() && cpts[i].cid == cpts[s].cid && cpts[i].face == cpts[s].face; ++i) {}
		// all cpts in [s,i) have the same {cid, face}
		crack_lines.resize(0);
		coll_obj const &cobj(coll_objects[cpts[s].cid]);
		cube_t const &cube(cobj);
		float const diameter(cube.get_bsphere_radius());
		
		for (unsigned j = s; j < i; ++j) { // generated cracks to the edge of the glass cube
			crack_point const &cpt1(cpts[j]);
			int const dim(cpt1.face >> 1), d1((dim+1)%3), d2((dim+2)%3);
			unsigned const ncracks(4); // one for each quadrant
			float const center(0.5*(cube.d[dim][0] + cube.d[dim][1]));
			float const x1(cpt1.pos[d1]), y1(cpt1.pos[d2]);
			rand_gen_t rgen;
			rgen.set_state(*(int *)&cpt1.orig_pos[d1], *(int *)&cpt1.orig_pos[d2]); // hash floats as ints	
			point epts[ncracks];

			for (unsigned n = 0; n < ncracks; ++n) {
				point epos;
				float min_dist_sq(0.0);

				for (unsigned attempt = 0; attempt < 4; ++attempt) {
					vector3d dir;
					dir[dim] = 0.0;
					dir[d1]  = rgen.rand_float()*((n&1) ? -1.0 : 1.0);
					dir[d2]  = rgen.rand_float()*((n&2) ? -1.0 : 1.0);
					point p1(cpt1.pos);
					p1[dim]  = center;
					point p2(p1 + dir.get_norm()*diameter);
					if (!do_line_clip(p1, p2, cube.d)) continue; // should never fail, and p1 should never change
					p2[dim]  = cpt1.pos[dim];

					for (vector<ray2d>::const_iterator c = crack_lines.begin(); c != crack_lines.end(); ++c) {
						float const x2(p2[d1]), x3(c->pts[0].x), x4(c->pts[1].x);
						if (max(x3, x4) < min(x1, x2) || max(x1, x2) < min(x3, x4)) continue;
						float const y2(p2[d2]), y3(c->pts[0].y), y4(c->pts[1].y);
						if (max(y3, y4) < min(y1, y2) || max(y1, y2) < min(y3, y4)) continue;
						float const denom((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1));
						if (fabs(denom) < TOLERANCE) continue;
						float const ub(((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3))/denom);
						if (ub < 0.0 || ub > 1.0)    continue;
						float const ua(((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3))/denom);
						if (ua < 0.0 || ua > 1.0)    continue;
						p2 = cpt1.pos + (p2 - cpt1.pos)*ua; // update intersection point
						if (attempt > 0 && p2p_dist_sq(cpt1.pos, p2) >= min_dist_sq) break;
					}
					float const dist_sq(p2p_dist_sq(cpt1.pos, p2));
					if (attempt == 0 || dist_sq < min_dist_sq) {epos = p2; min_dist_sq = dist_sq;}
				} // for attempt
				if (cpt1.pos != epos) {
					//beams.push_back(beam3d(0, NO_SOURCE, cpt1.pos, epos, cpt1.color, 0.05*cpt1.alpha));
					point pts[4] = {cpt1.pos, cpt1.pos, epos, epos};
					for (unsigned d = 0; d < 4; ++d) {pts[d][dim] = cube.d[dim][d == 1 || d == 2];}
					vector3d const normal(get_poly_norm(pts));
					colorRGBA color(cpt1.color);
					color.alpha *= min(0.5, 1.5*cobj.cp.color.alpha); // 1.5x glass alpha value
					qbd.add_quad_pts(pts, color, ((dot_product_ptv(normal, camera, pts[0]) < 0.0) ? -normal : normal));
				}
				epts[n] = epos;
			} // for n
			for (unsigned n = 0; n < ncracks; ++n) {
				crack_lines.push_back(ray2d(x1, y1, epts[n][d1], epts[n][d2]));
			}
		} // for j
	} // for i
}


void draw_cracks_and_decals() {

	if (!decals.any_active()) return;
	quad_batch_draw crack_qbd;
	create_and_draw_cracks(crack_qbd); // adds to beams
	map<int, quad_batch_draw> batches; // maps from {tid, is_black} to quad batches
	vector<pair<int, unsigned> > sorted_decals;
	line_tquad_draw_t blood_tqd;

	for (obj_vector_t<decal_obj>::const_iterator i = decals.begin(); i != decals.end(); ++i) {
		if (!i->status) continue;

		if (sphere_in_camera_view(i->get_pos(), i->radius, 0)) {
			sorted_decals.push_back(make_pair(-i->time, (i - decals.begin()))); // negate time, so largest time is first
		}
		i->maybe_draw_blood_trail(blood_tqd);
	}
	if (sorted_decals.empty() && crack_qbd.empty()) return;
	sort(sorted_decals.begin(), sorted_decals.end()); // sort by time, so that spraypaint works (later paint is drawn after/over earlier paint)

	for (unsigned i = 0; i < sorted_decals.size(); ++i) {
		decal_obj const &d(decals[sorted_decals[i].second]);
		d.draw(batches[(d.tid << 1) + (d.color == BLACK)]);
	}
	glDepthMask(GL_FALSE);
	enable_blend();
	shader_t black_shader, lighting_shader, bullet_shader;
	bool is_emissive(0);

	if (!blood_tqd.empty() || !crack_qbd.empty()) { // use normal lighting shader
		setup_smoke_shaders(lighting_shader, 0.01, 0, 1, 1, 1, 1, 1, 0, 1);
		lighting_shader.enable();
		select_texture(WHITE_TEX);
		blood_tqd.draw_tri_verts();
		blood_tqd.clear();
		crack_qbd.draw_and_clear();
	}
	for (map<int, quad_batch_draw>::const_iterator i = batches.begin(); i != batches.end(); ++i) {
		int const tid(i->first >> 1);
		bool const is_black(i->first & 1);

		if (tid == BULLET_D_TEX) { // use bullet shader
			if (!bullet_shader.is_setup()) {
				// see http://cowboyprogramming.com/2007/01/05/parallax-mapped-bullet-holes/
				bullet_shader.set_prefix("#define TEXTURE_ALPHA_MASK",  1); // FS
				bullet_shader.set_prefix("#define ENABLE_PARALLAX_MAP", 1); // FS
				setup_smoke_shaders(bullet_shader, 0.05, 0, 1, 1, 1, 1, 1, 0, 1, 1); // bump maps enabled
				bullet_shader.add_uniform_float("bump_tb_scale", -1.0); // invert the coordinate system (FIXME: something backwards?)
				bullet_shader.add_uniform_float("hole_depth", 0.2);
				bullet_shader.add_uniform_int("depth_map", 9);
				select_multitex(BULLET_N_TEX, 5, 0);
				select_multitex(BULLET_D_TEX, 9, 1);
			}
			bullet_shader.enable();
		}
		else if (is_black) { // use black shader
			if (!black_shader.is_setup()) {setup_smoke_shaders(black_shader, 0.01, 0, 1, 0, 0, 0, 1);} // no lighting
			black_shader.enable();
		}
		else { // use normal lighting shader
			if (!lighting_shader.is_setup()) {setup_smoke_shaders(lighting_shader, 0.01, 0, 1, 1, 1, 1, 1, 0, 1);}
			lighting_shader.enable();
			if (0 && !is_emissive) {lighting_shader.add_uniform_float("emissive_scale", 1.0); is_emissive = 1;} // make colors emissive
		}
		select_texture(tid);
		i->second.draw();
	} // for i
	disable_blend();
	glDepthMask(GL_TRUE);
	if (bullet_shader.is_setup()) {bullet_shader.enable(); bullet_shader.add_uniform_float("bump_tb_scale", 1.0);} // reset
	if (is_emissive) {lighting_shader.add_uniform_float("emissive_scale", 0.0);} // reset emissive
	black_shader.end_shader();
	lighting_shader.end_shader();
	bullet_shader.end_shader();
}


void setup_depth_trans_texture(shader_t &s, unsigned &depth_tid) {

	unsigned const depth_tu_id = 13;
	setup_depth_tex(s, depth_tu_id);
	depth_buffer_to_texture(depth_tid);
	bind_texture_tu(depth_tid, depth_tu_id);
	set_active_texture(0);
}

void draw_smoke_and_fires() {

	bool const use_depth_trans = 1;
	bool const have_part_clouds(part_clouds.any_active());
	if (!have_part_clouds && !fires.any_active()) return; // nothing to draw
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0.0, 0.0, use_depth_trans, 0, 0, 0); // no rain, snow, or reflections
	s.add_uniform_float("emissive_scale", 1.0); // make colors emissive
	set_multisample(0);
	unsigned depth_tid(0);
	if (use_depth_trans) {setup_depth_trans_texture(s, depth_tid);}

	if (have_part_clouds) { // smoke: slow when a lot of smoke is up close
		if (use_depth_trans) {s.add_uniform_float("depth_trans_bias", 0.1);}
		draw_part_clouds(part_clouds, SMOKE_PUFF_TEX);
	}
	order_vect_t fire_order;
	get_draw_order(fires, fire_order);
	
	if (!fire_order.empty()) {
		if (use_depth_trans) {s.add_uniform_float("depth_trans_bias", 0.01);}
		enable_blend();
		quad_batch_draw qbd;
		select_texture(FIRE_TEX);
		int last_in_smoke(-1);
		for (unsigned j = 0; j < fire_order.size(); ++j) {fires[fire_order[j].second].draw(qbd, last_in_smoke);}
		qbd.draw();
		set_std_blend_mode();
		disable_blend();
	}
	s.add_uniform_float("emissive_scale", 0.0); // reset
	s.end_shader();
	set_multisample(1);
	free_texture(depth_tid);
}


void add_camera_filter(colorRGBA const &color, unsigned time, int tid, unsigned ix, bool fades) {
	
	assert(ix < MAX_CFILTERS);
	if (color.alpha == 0.0) return;
	if (cfilters.size() <= ix) cfilters.resize(ix+1);
	cfilters[ix] = camera_filter(color, time, tid, fades);
}


void camera_filter::draw(bool apply_texture) {

	if (apply_texture) {select_texture(tid);} // use WHITE_TEX if tid < 0
	float const zval(-1.1*perspective_nclip), tan_val(tan(perspective_fovy/TO_DEG));
	float const y(0.5*zval*tan_val), x((y*window_width)/window_height);
	colorRGBA cur_color(color);
	if (fades) {cur_color.alpha *= float(time)/float(init_time);}
	cur_color.set_for_cur_shader();
	draw_tquad(x, y, zval);
}


void draw_camera_filters(vector<camera_filter> &cfs) {

	if (cfs.empty()) return;
	shader_t s;
	s.begin_simple_textured_shader();
	glDisable(GL_DEPTH_TEST);
	enable_blend();

	for (int i = (int)cfs.size()-1; i >= 0; --i) { // apply backwards
		if (cfs[i].time == 0) continue;
		cfs[i].draw();
		if ((int)cfs[i].time <= iticks) cfs[i].time = 0; else cfs[i].time -= iticks;
	}
	disable_blend();
	glEnable(GL_DEPTH_TEST);
	s.end_shader();
}


point world_space_to_screen_space(point const &pos) { // returns screen space normalized to [0.0, 1.0]

	double mats[2][16];
	fgGetMVM().get_as_doubles(mats[0]); // Model = MVM
	fgGetPJM().get_as_doubles(mats[1]); // Proj
	int const view[4] = {0, 0, 1, 1};
	vector3d_d pss;
	gluProject(pos.x, pos.y, pos.z, mats[0], mats[1], view, &pss.x, &pss.y, &pss.z);
	return point(pss);
}


void restore_prev_mvm_pjm_state() {

	fgMatrixMode(FG_PROJECTION);
	fgPopMatrix();
	fgMatrixMode(FG_MODELVIEW);
	fgPopMatrix();
}


bool is_sun_flare_visible() {

	if (!have_sun || light_factor < 0.4) return 0; // sun below the horizon, or doesn't exist
	point const cur_sun_pos(get_sun_pos());
	if (dot_product(cview_dir, (cur_sun_pos - get_camera_pos())) < 0.0) return 0; // sun behind the camera
	if (!sphere_in_camera_view(cur_sun_pos, 4.0*sun_radius, 2)) return 0; // use larger radius to include the flare/halo
	return 1;
}


float const spark_t::radius = 0.0;


void spark_t::draw(quad_batch_draw &qbd) const {

	point const camera(get_camera_pos());
	qbd.add_billboard((pos + (camera - pos).get_norm()*0.02), camera, up_vector, c, 0.8*s, 0.8*s);
}


void draw_sparks() {

	if (sparks.empty()) return;
	ensure_filled_polygons();
	enable_blend();
	set_additive_blend_mode();
	shader_t s;
	s.begin_simple_textured_shader(0.01);
	select_texture(FLARE2_TEX);
	quad_batch_draw qbd;
	draw_objects(sparks, qbd);
	qbd.draw();
	s.end_shader();
	set_std_blend_mode();
	disable_blend();
	reset_fill_mode();
	sparks.clear();
}


void draw_projectile_effects() {

	draw_blasts();
	draw_beams();
	draw_sparks();
	explosion_part_man.draw(0.002, -1); // untextured
	water_part_man.draw(); // not really a projectile effect, but it's drawn with them
}


struct splash_ring_t {

	point pos;
	float size;
	colorRGBA color;

	splash_ring_t(point const &pos_, float size_, colorRGBA const &color_) : pos(pos_), size(size_), color(color_) {}

	void draw(shader_t &shader) const {
		unsigned const num_rings(min(10U, (unsigned)ceil(size)));
		float radius(min(size, 0.025f));
		float const dr(0.5*radius);
		unsigned const ndiv(max(3, min(N_CYL_SIDES, int(1000.0*radius/max(TOLERANCE, distance_to_camera(pos))))));
		shader.set_cur_color(color);

		for (unsigned i = 0; i < num_rings; ++i) {
			draw_circle_normal((radius - 0.5*dr), radius, ndiv, 0, pos);
			radius += dr;
		}
	}
};


vector<splash_ring_t> splashes;


void draw_splash(float x, float y, float z, float size, colorRGBA color) { // queue it up for drawing

	assert(size >= 0.0);
	if (DISABLE_WATER || !(display_mode & 0x04)) return;
	if (size == 0.0 || temperature <= W_FREEZE_POINT) return;
	if (size > 0.1) size = sqrt(10.0*size)/10.0;
	select_liquid_color(color, get_xpos(x), get_ypos(y));
	splashes.push_back(splash_ring_t(point(x, y, z+0.001), size, color));
}


void draw_splashes() {

	if (splashes.empty()) return;
	shader_t s;
	s.begin_untextured_lit_glcolor_shader();
	for (auto i = splashes.begin(); i != splashes.end(); ++i) {i->draw(s);}
	s.end_shader();
	splashes.clear(); // only last one frame
}


void draw_framerate(float val) {

	char text[32];
	sprintf(text, "%i", round_fp(val));
	float const ar(((float)window_width)/((float)window_height));
	draw_text(WHITE, -0.011*ar, -0.011, -2.0*NEAR_CLIP, text);
}


void draw_compass_and_alt() { // and temperature

	char text[64];
	float const aspect_ratio((float)window_width/(float)window_height);
	string const dirs[8] = {"N", "NW", "W", "SW", "S", "SE", "E", "NE"};
	sprintf(text, "Loc: (%3.2f, %3.2f, %3.2f)", (camera_origin.x+(xoff2-xoff)*DX_VAL), (camera_origin.y+(yoff2-yoff)*DY_VAL), camera_origin.z);
	draw_text(YELLOW, -0.005*aspect_ratio, -0.01, -0.02, text);
	float const theta(safe_acosf(cview_dir.x)*TO_DEG);
	int const octant(int(((cview_dir.y < 0) ? (360.0 - theta) : theta)/45.0 + 22.5)&7);
	sprintf(text, "%s", dirs[octant].c_str());
	draw_text(YELLOW, 0.005*aspect_ratio, -0.01, -0.02, text);
	sprintf(text, "Temp: %iC", int(temperature));
	draw_text(YELLOW, 0.007*aspect_ratio, -0.01, -0.02, text);
}

