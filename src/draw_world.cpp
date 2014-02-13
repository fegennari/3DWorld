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


// Global Variables
float sun_radius, moon_radius, earth_radius, brightness(1.0);
colorRGBA cur_ambient(BLACK), cur_diffuse(BLACK);
point sun_pos, moon_pos;
point gl_light_positions[8] = {all_zeros};
point const earth_pos(-15.0, -8.0, 21.0);
sky_pos_orient cur_spo(point(0,0,0),1,0,0);
vector3d up_norm(plus_z);
vector<camera_filter> cfilters;
pt_line_drawer bubble_pld;

extern bool have_sun, using_lightmap, has_dl_sources, has_spotlights, has_line_lights, smoke_exists, two_sided_lighting;
extern bool group_back_face_cull, have_indir_smoke_tex, combined_gu;
extern int is_cloudy, iticks, frame_counter, display_mode, show_fog, num_groups, xoff, yoff;
extern int window_width, window_height, game_mode, enable_fsource, draw_model, camera_mode, DISABLE_WATER;
extern unsigned smoke_tid, dl_tid, num_stars, create_voxel_landscape;
extern float zmin, light_factor, fticks, perspective_fovy, perspective_nclip, cobj_z_bias;
extern float temperature, atmosphere, zbottom, indir_vert_offset;
extern point light_pos, mesh_origin, flow_source, surface_pos;
extern vector3d wind;
extern colorRGB const_indir_color, ambient_lighting_scale;
extern colorRGBA bkg_color, sun_color;
extern vector<spark_t> sparks;
extern vector<star> stars;
extern vector<beam3d> beams;
extern obj_group obj_groups[];
extern coll_obj_group coll_objects;
extern obj_type object_types[];
extern obj_vector_t<bubble> bubbles;
extern obj_vector_t<particle_cloud> part_clouds;
extern cloud_manager_t cloud_manager;
extern obj_vector_t<fire> fires;
extern obj_vector_t<decal_obj> decals;
extern water_particle_manager water_part_man;
extern cube_t cur_smoke_bb;
extern vector<portal> portals;
extern vector<obj_draw_group> obj_draw_groups;



void set_fill_mode() {
	glPolygonMode(GL_FRONT_AND_BACK, ((draw_model == 0) ? GL_FILL : GL_LINE));
}

int get_universe_ambient_light() {
	return ((world_mode == WMODE_UNIVERSE) ? GL_LIGHT1 : GL_LIGHT3);
}


void set_colors_and_enable_light(int light, float const ambient[4], float const diffuse[4]) {

	assert(light >= GL_LIGHT0 && light <= GL_LIGHT7);
	glEnable(light);
	glLightfv(light, GL_AMBIENT, ambient);
	glLightfv(light, GL_DIFFUSE, diffuse);
}


void clear_colors_and_disable_light(int light) {

	assert(light >= GL_LIGHT0 && light <= GL_LIGHT7);
	float const ad[4] = {0.0, 0.0, 0.0, 0.0};
	glDisable(light);
	glLightfv(light, GL_AMBIENT, ad);
	glLightfv(light, GL_DIFFUSE, ad);
}


void set_gl_light_pos(int light, point const &pos, float w) {

	assert(light >= GL_LIGHT0 && light <= GL_LIGHT7);
	float const position[4] = {pos.x, pos.y, pos.z, w};
	glLightfv(light, GL_POSITION, position);
	gl_light_positions[light - GL_LIGHT0] = pos;
}


void set_color_alpha(colorRGBA color, float alpha) {

	color.alpha *= alpha;
	colorRGBA(0.0, 0.0, 0.0, color.alpha).do_glColor(); // sets alpha component
	set_color_a(BLACK);
	// FIXME: replace black with near black to fix a color setting bug that sometimes leaves the object as white when black is specified
	set_color_d((color == BLACK) ? colorRGBA(0.001, 0.001, 0.001, color.alpha) : color);
}


void draw_camera_weapon(bool want_has_trans) {

	if (!game_mode || weap_has_transparent(CAMERA_ID) != want_has_trans) return;
	shader_t s;
	setup_smoke_shaders(s, 0.0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1);
	draw_weapon_in_hand(-1);
	s.end_shader();
}


void set_specular(float specularity, float shininess) {

	static float last_shiny(-1.0), last_spec(-1.0);
	if (is_cloudy && world_mode != WMODE_UNIVERSE) specularity *= 0.5;

	if (specularity != last_spec) { // This materialfv stuff seems to take some time, so only set if changed since last call
		float mat_specular[]  = {specularity, specularity, specularity, 1.0};
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
		last_spec = specularity;
	}
	if (shininess != last_shiny) {
		float mat_shininess[] = {max(0.0f, min(128.0f, shininess))};
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		last_shiny = shininess;
    }
}


void calc_cur_ambient_diffuse() {

	float a[4], d[4], lval[4];
	unsigned ncomp(0);
	cur_ambient = cur_diffuse = BLACK;

	for (unsigned i = 0; i < 8; ++i) { // max of 8 lights (GL_LIGHT0 - GL_LIGHT7): sun, moon, lightning
		int const light(GL_LIGHT0 + i); // should be sequential

		if (glIsEnabled(light)) {
			float atten(1.0);
			glGetLightfv(light, GL_AMBIENT, a);
			glGetLightfv(light, GL_DIFFUSE, d);
			glGetLightfv(light, GL_POSITION, lval);
			if (lval[3] != 0.0) glGetLightfv(light, GL_CONSTANT_ATTENUATION, &atten); // point light source only
			assert(atten > 0.0);
			UNROLL_3X(cur_ambient[i_] += a[i_]/atten; cur_diffuse[i_] += d[i_]/atten;)
			//cout << "A: "; cur_ambient.print(); cout << "  D: "; cur_diffuse.print(); cout << endl;
			++ncomp;
		}
	}
	if (ncomp > 0) {
		float const cscale(0.5 + 0.5/ncomp);
		cur_ambient       = cur_ambient.modulate_with(ambient_lighting_scale);
		cur_ambient      *= cscale; // only really valid for sun and moon
		cur_diffuse      *= cscale;
		cur_ambient.alpha = 1.0;
		cur_diffuse.alpha = 1.0;
	}
}


void upload_mvm_to_shader(shader_t &s, char const *const var_name) {

	float mvm[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mvm);
	s.add_uniform_matrid_4x4(var_name, mvm, 0);
}


void set_dlights_booleans(shader_t &s, bool enable, int shader_type) {

	if (!enable)         {s.set_prefix("#define NO_DYNAMIC_LIGHTS", shader_type);} // if we're not even enabling dlights
	if (has_spotlights)  {s.set_prefix("#define HAS_SPOTLIGHTS",    shader_type);}
	if (has_line_lights) {s.set_prefix("#define HAS_LINE_LIGHTS",   shader_type);}
	s.set_bool_prefix("enable_dlights", (enable && dl_tid > 0 && has_dl_sources), shader_type);
}


void common_shader_block_pre(shader_t &s, bool &dlights, bool &use_shadow_map, bool &indir_lighting, float min_alpha) {

	use_shadow_map &= shadow_map_enabled();
	indir_lighting &= have_indir_smoke_tex;
	dlights        &= (dl_tid > 0 && has_dl_sources);
	s.check_for_fog_disabled();
	if (min_alpha == 0.0) {s.set_prefix("#define NO_ALPHA_TEST", 1);} // FS
	s.set_bool_prefixes("indir_lighting", indir_lighting, 3); // VS/FS
	s.set_bool_prefix("use_shadow_map", use_shadow_map, 1); // FS
	set_dlights_booleans(s, dlights, 1); // FS
}


void set_indir_lighting_block(shader_t &s, bool use_smoke_indir) {

	if (use_smoke_indir && smoke_tid) {set_3d_texture_as_current(smoke_tid, 1);}
	s.add_uniform_int("smoke_and_indir_tex", 1);
	s.add_uniform_float("half_dxy", HALF_DXY);
	s.add_uniform_float("indir_vert_offset", indir_vert_offset);
	colorRGB const indir_color((have_indir_smoke_tex && world_mode == WMODE_GROUND) ? colorRGB(0.0, 0.0, 0.0) : const_indir_color);
	s.add_uniform_color("const_indir_color", indir_color);
}


void common_shader_block_post(shader_t &s, bool dlights, bool use_shadow_map, bool use_smoke_indir, float min_alpha) {

	s.setup_scene_bounds();
	s.setup_fog_scale(); // fog scale for the case where smoke is disabled
	if (dlights) setup_dlight_textures(s);
	set_indir_lighting_block(s, use_smoke_indir);
	s.add_uniform_int("tex0", 0);
	s.add_uniform_float("min_alpha", min_alpha);
	if (use_shadow_map) set_smap_shader_for_all_lights(s, cobj_z_bias);
	set_active_texture(0);
}


void set_smoke_shader_prefixes(shader_t &s, int use_texgen, bool keep_alpha, bool direct_lighting,
	bool smoke_enabled, bool has_lt_atten, int use_bmap, bool use_spec_map, bool use_mvm, bool use_tsl)
{
	s.set_int_prefix ("use_texgen",      use_texgen,      0); // VS
	s.set_bool_prefix("keep_alpha",      keep_alpha,      1); // FS
	s.set_bool_prefix("direct_lighting", direct_lighting, 1); // FS
	s.set_bool_prefix("do_lt_atten",     has_lt_atten,    1); // FS
	s.set_bool_prefix("two_sided_lighting",  use_tsl,     1); // FS
	s.set_bool_prefix("use_world_space_mvm", use_mvm,     0); // VS
	if (use_spec_map) {s.set_prefix("#define USE_SPEC_MAP", 1);} // FS

	if (smoke_enabled) {
		// Note: dynamic_smoke_shadows applies to light0 only
		// Note: dynamic_smoke_shadows still uses the visible smoke bbox, so if you can't see smoke it won't cast a shadow
		for (unsigned d = 0; d < 2; ++d) { // VS/FS
			if (DYNAMIC_SMOKE_SHADOWS) {s.set_prefix("#define DYNAMIC_SMOKE_SHADOWS", d);}
			s.set_prefix("#define SMOKE_ENABLED", d);
		}
	}
	for (unsigned i = 0; i < 2; ++i) {
		if (use_bmap     ) {s.set_prefix("#define USE_BUMP_MAP",       i);} // VS/FS
		if (use_bmap == 2) {s.set_prefix("#define USE_TANGENT_VECTOR", i);} // VS/FS
	}
	s.setup_enabled_lights(8, 2); // FS
}


// texture units used: 0: object texture, 1: smoke/indir lighting texture, 2-4 dynamic lighting, 5: bump map, 6-7 shadow map, 8: specular map, 9: depth map, 10: burn mask
// use_texgen: 0 = use texture coords, 1 = use standard texture gen matrix, 2 = use custom shader tex0_s/tex0_t, 3 = use vertex id for texture
void setup_smoke_shaders(shader_t &s, float min_alpha, int use_texgen, bool keep_alpha, bool indir_lighting, bool direct_lighting, bool dlights, bool smoke_en,
	bool has_lt_atten, bool use_smap, int use_bmap, bool use_spec_map, bool use_mvm, bool force_tsl, bool use_light_colors, float burn_offset)
{
	bool const use_burn_mask(burn_offset > -1.0);
	smoke_en &= (have_indir_smoke_tex && smoke_exists && smoke_tid > 0);
	if (use_light_colors) {s.set_prefix("#define USE_LIGHT_COLORS", 1);} // FS
	if (use_burn_mask   ) {s.set_prefix("#define APPLY_BURN_MASK",  1);} // FS
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, min_alpha);
	set_smoke_shader_prefixes(s, use_texgen, keep_alpha, direct_lighting, smoke_en, has_lt_atten, use_bmap, use_spec_map, use_mvm, force_tsl);
	s.set_vert_shader("texture_gen.part+line_clip.part*+bump_map.part+indir_lighting.part+no_lt_texgen_smoke");
	s.set_frag_shader("fresnel.part*+linear_fog.part+bump_map.part+spec_map.part+ads_lighting.part*+dynamic_lighting.part*+shadow_map.part*+line_clip.part*+indir_lighting.part+black_body_burn.part+textured_with_smoke");
	s.begin_shader();

	if (use_texgen == 2) {
		s.register_attrib_name("tex0_s", TEX0_S_ATTR);
		s.register_attrib_name("tex0_t", TEX0_T_ATTR);
	}
	if (use_bmap)     s.add_uniform_int("bump_map", 5);
	if (use_spec_map) s.add_uniform_int("spec_map", 8);
	s.add_uniform_float("base_color_scale", (use_light_colors ? 0.0 : 1.0)); // hack to force usage of material properties instead of color
	common_shader_block_post(s, dlights, use_smap, (smoke_en || indir_lighting), min_alpha);
	float const step_delta_scale(get_smoke_at_pos(get_camera_pos()) ? 1.0 : 2.0);
	s.add_uniform_float_array("smoke_bb", &cur_smoke_bb.d[0][0], 6);
	s.add_uniform_float("step_delta", step_delta_scale*HALF_DXY);
	if (use_mvm ) {upload_mvm_to_shader(s, "world_space_mvm");}
	if (smoke_en) {s.add_uniform_color("smoke_color", colorRGB(GRAY));}

	if (use_burn_mask) {
		s.add_uniform_float("burn_tex_scale", 0.05); // FIXME: hard-coded
		s.add_uniform_float("burn_offset", burn_offset);
		s.add_uniform_int("burn_mask", 10);
		select_multitex(DISINT_TEX, 10); // PLASMA_TEX?
	}
}


void set_tree_branch_shader(shader_t &s, bool direct_lighting, bool dlights, bool use_smap) {

	bool indir_lighting(0);
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, 0.0);
	set_smoke_shader_prefixes(s, 0, 0, direct_lighting, 0, 0, 0, 0, 0, 0);
	s.set_vert_shader("texture_gen.part+line_clip.part*+bump_map.part+indir_lighting.part+no_lt_texgen_smoke");
	s.set_frag_shader("fresnel.part*+linear_fog.part+bump_map.part+ads_lighting.part*+dynamic_lighting.part*+shadow_map.part*+line_clip.part*+indir_lighting.part+textured_with_smoke");
	s.begin_shader();
	common_shader_block_post(s, dlights, use_smap, 0, 0.0);
	check_gl_error(400);
}


// texture units used: 0,8,15: object texture, 1: indir lighting texture, 2-4: dynamic lighting, 5: 3D noise texture, 6-7: shadow map, 9-14: tree leaf textures | 9: AO texture, 10: voxel shadow texture
void setup_procedural_shaders(shader_t &s, float min_alpha, bool indir_lighting, bool dlights, bool use_smap,
	bool use_noise_tex, bool z_top_test, float tex_scale, float noise_scale, float tex_mix_saturate)
{
	common_shader_block_pre(s, dlights, use_smap, indir_lighting, min_alpha);
	s.set_bool_prefix("use_noise_tex",  use_noise_tex,  1); // FS
	s.set_bool_prefix("z_top_test",     z_top_test,     1); // FS
	s.setup_enabled_lights(2, 2); // FS; only 2, but could be up to 8 later
	s.set_vert_shader("indir_lighting.part+procedural_gen");
	s.set_frag_shader("linear_fog.part+ads_lighting.part*+dynamic_lighting.part*+shadow_map.part*+triplanar_texture.part+procedural_texture.part+indir_lighting.part+voxel_texture.part+procedural_gen");
	s.begin_shader();
	common_shader_block_post(s, dlights, use_smap, indir_lighting, min_alpha);
	s.add_uniform_int("tex1",    8);
	s.add_uniform_int("tex_top", 15); // not used in all cases
	s.add_uniform_float("tex_scale", tex_scale);

	if (use_noise_tex) {
		s.add_uniform_int("noise_tex", 5); // does this need an enable option?
		s.add_uniform_float("noise_scale", noise_scale);
		s.add_uniform_float("tex_mix_saturate", tex_mix_saturate);
	}
}


void setup_object_render_data() {

	RESET_TIME;
	bool const TIMETEST(0);
	calc_cur_ambient_diffuse();
	distribute_smoke();
	if (TIMETEST) {PRINT_TIME("1 Distribute Smoke");}
	upload_smoke_indir_texture();
	if (TIMETEST) {PRINT_TIME("2 Upload Smoke");}
	add_dynamic_lights_ground();
	if (TIMETEST) {PRINT_TIME("3 Add Dlights");}
	cube_t const dlight_bounds(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, get_zval_min(), get_zval_max());
	upload_dlights_textures(dlight_bounds); // get_scene_bounds()
	if (TIMETEST) {PRINT_TIME("4 Dlights Textures");}
	get_occluders();
	if (TIMETEST) {PRINT_TIME("5 Get Occluders");}
}


void end_group(int &last_group_id) {

	if (last_group_id < 0) return;
	assert((unsigned)last_group_id < obj_draw_groups.size());
	obj_draw_groups[last_group_id].end_render();
	if (group_back_face_cull) glDisable(GL_CULL_FACE);
	last_group_id = -1;
}


// should always have draw_solid enabled on the first call for each frame
void draw_coll_surfaces(bool draw_solid, bool draw_trans) {

	//RESET_TIME;
	assert(draw_solid || draw_trans);
	static vector<pair<float, int> > draw_last;
	if (coll_objects.empty() || coll_objects.drawn_ids.empty() || world_mode != WMODE_GROUND) return;
	if (!draw_solid && draw_last.empty() && (!smoke_exists || portals.empty())) return; // nothing transparent to draw
	set_fill_mode();
	// Note: in draw_solid mode, we could call get_shadow_triangle_verts() on occluders to do a depth pre-pass here, but that doesn't seem to be more efficient
	glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
	glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
	set_color_a(BLACK);
	set_specular(0.0, 1.0);
	bool has_lt_atten(draw_trans && !draw_solid && coll_objects.has_lt_atten);
	// Note: enable direct_lighting if processing sun/moon shadows here
	float burn_offset(-1.0);
	//if (display_mode & 0x10) {burn_offset = 0.002*(frame_counter%1000) - 1.0;}
	shader_t s;
	setup_smoke_shaders(s, 0.0, 2, 0, 1, 1, 1, 1, has_lt_atten, 1, 0, 0, 0, two_sided_lighting, 0, burn_offset);
	if (!s.is_setup()) has_lt_atten = 0; // shaders disabled
	int last_tid(-1), last_group_id(-1);
	vector<vert_wrap_t> portal_verts;
	vector<vert_norm> poly_verts;
	
	if (draw_solid) {
		vector<pair<float, int> > large_cobjs;
		draw_last.resize(0);

		for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
			unsigned cix(*i);
			assert(cix < coll_objects.size());
			coll_obj const &c(coll_objects[cix]);
			assert(c.cp.draw);
			if (c.no_draw()) continue; // can still get here sometimes

			if (c.is_big_occluder() && c.group_id < 0) {
				float const dist(distance_to_camera(c.get_center_pt()));

				if (c.get_area() > 0.01*dist*dist) { // increases CPU time but decreases GPU time
					if (camera_pdu.cube_visible(c)) {large_cobjs.push_back(make_pair(dist, *i));}
					continue;
				}
			}
			if (c.is_semi_trans()) { // slow when polygons are grouped
				float dist(distance_to_camera(c.get_center_pt()));

				if (c.type == COLL_SPHERE) { // distance to surface closest to the camera
					dist -= c.radius;
				}
				else if (c.type == COLL_CYLINDER || c.type == COLL_CYLINDER_ROT) { // approx distance to surface closest to the camera
					dist -= min(0.5*(c.radius + c.radius2), 0.5*p2p_dist(c.points[0], c.points[1]));
				}
				draw_last.push_back(make_pair(-dist, cix)); // negative distance
			}
			else {
				c.draw_cobj(cix, last_tid, last_group_id, poly_verts, s); // i may not be valid after this call
				
				if (cix != *i) {
					assert(cix > *i);
					i = std::lower_bound(i, coll_objects.drawn_ids.end(), cix);
				}
			}
		} // for i
		end_group(last_group_id);
		sort(large_cobjs.begin(), large_cobjs.end()); // sort front to back for early z culling

		for (vector<pair<float, int> >::const_iterator i = large_cobjs.begin(); i != large_cobjs.end(); ++i) {
			unsigned cix(i->second);
			coll_objects[cix].draw_cobj(cix, last_tid, last_group_id, poly_verts, s);
		}
	} // end draw solid
	if (draw_trans) { // called second
		if (smoke_exists) {
			for (unsigned i = 0; i < portals.size(); ++i) {
				if (!portals[i].is_visible()) continue;
				draw_last.push_back(make_pair(-distance_to_camera(portals[i].get_center_pt()), -(int)(i+1)));
			}
		}
		sort(draw_last.begin(), draw_last.end()); // sort back to front for alpha blending
		enable_blend();
		int ulocs[3] = {0};
		float last_light_atten(-1.0), last_refract_ix(0.0); // set to invalid values to start
		bool in_portal(0);

		if (has_lt_atten) {
			ulocs[0] = s.get_uniform_loc("light_atten");
			ulocs[1] = s.get_uniform_loc("cube_bb"    );
			ulocs[2] = s.get_uniform_loc("refract_ix" );
			assert(ulocs[0] && ulocs[1] && ulocs[2]);
		}
		for (unsigned i = 0; i < draw_last.size(); ++i) {
			int const ix(draw_last[i].second);

			if (ix < 0) { // portal
				end_group(last_group_id);

				if (has_lt_atten && last_light_atten != 0.0) {
					s.set_uniform_float(ulocs[0], 0.0);
					last_light_atten = 0.0;
				}
				if (has_lt_atten && last_refract_ix != 1.0) {
					s.set_uniform_float(ulocs[2], 1.0);
					last_refract_ix = 1.0;
				}
				if (!in_portal) {portal::pre_draw(portal_verts); in_portal = 1;}
				unsigned const pix(-(ix+1));
				assert(pix < portals.size());
				portals[pix].draw(portal_verts);
			}
			else { // cobj
				if (in_portal) {portal::post_draw(portal_verts); in_portal = 0;}
				unsigned cix(ix);
				assert(cix < coll_objects.size());
				coll_obj const &c(coll_objects[cix]);
				
				if (has_lt_atten) { // we only support cubes for now (Note: may not be compatible with groups)
					float const light_atten((c.type == COLL_CUBE) ? c.cp.light_atten : 0.0);

					if (light_atten != last_light_atten) {
						s.set_uniform_float(ulocs[0], light_atten);
						last_light_atten = light_atten;
					}
					if (c.cp.refract_ix != last_refract_ix) {
						s.set_uniform_float(ulocs[2], c.cp.refract_ix);
						last_refract_ix = c.cp.refract_ix;
					}
					if (light_atten > 0.0) s.set_uniform_float_array(ulocs[1], (float const *)c.d, 6);
				}
				c.draw_cobj(cix, last_tid, last_group_id, poly_verts, s);
				assert(cix == ix); // should not have changed
			}
		} // for i
		if (in_portal) {portal::post_draw(portal_verts);}
		end_group(last_group_id);
		disable_blend();
		draw_last.resize(0);
	} // end draw_trans
	s.end_shader();
	set_specular(0.0, 1.0);
	//if (draw_solid) PRINT_TIME("Final Draw");
}


bool portal::is_visible() const {

	point center;
	float rad;
	polygon_bounding_sphere(pts, 4, 0.0, center, rad);
	if (normal != zero_vector && dot_product_ptv(normal, get_camera_pos(), center) < 0.0) return 0; // back facing
	return sphere_in_camera_view(center, rad, 2);
}

void portal::pre_draw(vector<vert_wrap_t> &verts) {

	float const scale[2] = {0.0, 0.0}, xlate[2] = {0.0, 0.0};
	select_texture(WHITE_TEX, 0);
	setup_polygon_texgen(plus_z, scale, xlate, zero_vector); // doesn't matter as long as it's set to something
	ALPHA0.do_glColor();
	assert(verts.empty());
}

void portal::post_draw(vector<vert_wrap_t> &verts) {

	draw_verts(verts, GL_QUADS);
	verts.clear();
}

void portal::draw(vector<vert_wrap_t> &verts) const {

	for (unsigned i = 0; i < 4; ++i) {verts.push_back(pts[i]);}
};


void draw_stars(float alpha) {

	assert(num_stars <= stars.size());
	if (alpha <= 0.0) return;
	colorRGBA color(BLACK), bkg;
	UNROLL_3X(bkg[i_] = (1.0 - alpha)*bkg_color[i_];)
	glPushMatrix();
	if (camera_mode == 1) {translate_to(surface_pos);}
	set_color(BLACK);
	enable_blend();
	glPointSize(2.0);
	glDisable(GL_DEPTH_TEST);
	shader_t s;
	s.begin_color_only_shader();
	vector<vert_color> pts;
	pts.reserve(num_stars);

	for (unsigned i = 0; i < num_stars; ++i) {
		if ((rand()%400) == 0) continue; // flicker out

		for (unsigned j = 0; j < 3; ++j) {
			float const c(stars[i].color[j]*stars[i].intensity);
			color[j] = ((alpha >= 1.0) ? c : (alpha*c + bkg[j]));
		}
		pts.push_back(vert_color(stars[i].pos, color));
	}
	draw_verts(pts, GL_POINTS);
	s.end_shader();
	glEnable(GL_DEPTH_TEST);
	glPointSize(1.0);
	disable_blend();
	glPopMatrix();
}


void draw_sun() {

	point const pos(get_sun_pos());
	if (!have_sun || !sphere_in_camera_view(pos, sun_radius, 1)) return;
	//select_texture(SUN_TEX);
	glDisable(GL_LIGHTING);
	colorRGBA color(SUN_C);
	apply_red_sky(color);
	color.do_glColor();
	draw_subdiv_sphere(pos, sun_radius, N_SPHERE_DIV, 1, 0);
	glEnable(GL_LIGHTING);
	//glDisable(GL_TEXTURE_2D);
}


void draw_moon() {

	if (world_mode == WMODE_GROUND && show_fog) return; // don't draw when there is fog
	point const pos(get_moon_pos());
	if (!sphere_in_camera_view(pos, moon_radius, 1)) return;
	set_color(WHITE);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);

	if (have_sun) {
		float const ambient[4] = {0.05, 0.05, 0.05, 1.0}, diffuse[4] = {1.0, 1.0, 1.0, 1.0};
		set_gl_light_pos(GL_LIGHT4, get_sun_pos(), 0.0);
		set_colors_and_enable_light(GL_LIGHT4, ambient, diffuse);
	}
	select_texture(MOON_TEX);
	draw_subdiv_sphere(pos, moon_radius, N_SPHERE_DIV, 1, 0);
	glDisable(GL_TEXTURE_2D);
	if (light_factor < 0.6) glEnable(GL_LIGHT1); // moon
	if (light_factor > 0.4) glEnable(GL_LIGHT0); // sun
	glDisable(GL_LIGHT4);

	if (light_factor >= 0.4) { // fade moon into background color when the sun comes up
		colorRGBA color(bkg_color);
		color.alpha = 5.0*(light_factor - 0.4);
		glDisable(GL_LIGHTING);
		enable_blend();
		color.do_glColor();
		draw_subdiv_sphere(pos, 1.1*moon_radius, N_SPHERE_DIV, 0, 0);
		glEnable(GL_LIGHTING);
		disable_blend();
	}
}


// for some reason the texture is backwards, so we mirrored the image of the earth
void draw_earth() {

	if (show_fog) return; // don't draw when there is fog
	point pos(mesh_origin + earth_pos);
	if (camera_mode == 1) pos += surface_pos;
	static float rot_angle(0.0);

	if (sphere_in_camera_view(pos, earth_radius, 1)) {
		set_fill_mode();
		select_texture(EARTH_TEX);
		set_color(WHITE);
		glPushMatrix();
		translate_to(pos);
		glRotatef(67.0, 0.6, 0.8, 0.0);
		glRotatef(rot_angle, 0.0, 0.0, 1.0);
		glRotatef(180.0, 1.0, 0.0, 0.0);
		draw_sphere_vbo(all_zeros, earth_radius, N_SPHERE_DIV, 1);
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);
	}
	rot_angle += 0.2*fticks;
}


void draw_stationary_earth(float radius) {

	set_fill_mode();
	select_texture(EARTH_TEX);
	set_color(WHITE);
	draw_subdiv_sphere(all_zeros, radius, N_SPHERE_DIV, 1, 0);
	glDisable(GL_TEXTURE_2D);
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

	colorRGBA color(brightness, brightness, brightness, atmosphere);
	apply_red_sky(color);
	return color;
}


void get_avg_sky_color(colorRGBA &avg_color) {

	colorRGBA cloud_color(get_cloud_color());
	cloud_color.alpha = 1.0;
	blend_color(avg_color, cloud_color, bkg_color, 0.5, 1);
}


float get_cloud_density(point const &pt, vector3d const &dir) { // optimize?

	if (atmosphere == 0.0) return 0.0;
	point lsint;
	if (!line_sphere_int(-dir, pt, cur_spo.center, cur_spo.radius, lsint, 0)) return 0.0; // shouldn't get here?
	vector3d const vdir(lsint - cur_spo.center);
	return atmosphere*get_texture_component(CLOUD_TEX, (vdir.x*cur_spo.radius_inv + cur_spo.dx), (vdir.y*cur_spo.radius_inv + cur_spo.dy), 3); // cloud alpha
}


void draw_sky(int order) { // FIXME SHADERS: uses fixed function pipeline

	if (atmosphere < 0.01) return; // no atmosphere
	set_specular(0.0, 1.0);
	float radius(0.55*(FAR_CLIP+X_SCENE_SIZE));
	point center((camera_mode == 1) ? surface_pos : mesh_origin);
	center.z -= 0.727*radius;
	if ((distance_to_camera(center) > radius) != order) return;
	colorRGBA const cloud_color(get_cloud_color());

	static float sky_rot_xy[2] = {0.0, 0.0}; // x, y
	float const wmag(sqrt(wind.x*wind.x + wind.y*wind.y));

	if (wmag > TOLERANCE) {
		for (unsigned d = 0; d < 2; ++d) {
			sky_rot_xy[d] -= fticks*CLOUD_WIND_SPEED*(wmag + 0.5*WIND_ADJUST)*wind[d]/wmag;
		}
	}
	cur_spo = sky_pos_orient(center, radius, sky_rot_xy[0], sky_rot_xy[1]);
	int const light(GL_LIGHT4);
	set_fill_mode();
	enable_blend();

	if (have_sun && light_factor > 0.4) { // sun lighting of clouds
		float diffuse[4], ambient[4];
		point lpos(get_sun_pos()), lsint;
		vector3d const sun_v((get_camera_pos() - lpos).get_norm());
		if (line_sphere_int(sun_v, lpos, center, radius, lsint, 1)) lpos = lsint;
		
		for (unsigned i = 0; i < 4; ++i) { // even alpha?
			diffuse[i] = 1.0*sun_color[i];
			ambient[i] = 0.5*sun_color[i];
		}
		set_gl_light_pos(light, lpos, 1.0); // w - point light source
		set_colors_and_enable_light(light, ambient, diffuse);
		glLightf(light, GL_CONSTANT_ATTENUATION,  0.0);
		glLightf(light, GL_LINEAR_ATTENUATION,    0.01);
		glLightf(light, GL_QUADRATIC_ATTENUATION, 0.01);
	}
	if (have_sun && light_factor > 0.4) { // draw horizon
		glDisable(GL_LIGHTING);
		colorRGBA horizon_color;
		float const blend_val(atmosphere*CLIP_TO_01(10.0f*(light_factor - 0.4f)));
		blend_color(horizon_color, WHITE, ALPHA0, blend_val, 1);
		horizon_color.alpha *= 0.5;
		apply_red_sky(horizon_color);
		horizon_color.do_glColor();
		select_texture(GRADIENT_TEX);
		draw_sphere_vbo(center, 1.05*radius, N_SPHERE_DIV, 1);
		glEnable(GL_LIGHTING);
	}
	select_texture(CLOUD_TEX);

	// change S and T parameters to map sky texture into the x/y plane with translation based on wind/rot
	bool const depth_clamp_enabled(glIsEnabled(GL_DEPTH_CLAMP) != 0);
	if (depth_clamp_enabled) {glDisable(GL_DEPTH_CLAMP);}
	glEnable(GL_TEXTURE_GEN_S); glEnable(GL_TEXTURE_GEN_T);
	setup_texgen(1.0/radius, 1.0/radius, (sky_rot_xy[0] - center.x/radius), (sky_rot_xy[1] - center.y/radius)); // GL_EYE_LINEAR
	set_color_a(cloud_color);
	set_color_d(cloud_color); // disable lighting (BLACK)?
	draw_subdiv_sphere(center, radius, (3*N_SPHERE_DIV)/2, zero_vector, NULL, 0, 1);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_TEXTURE_GEN_S); glDisable(GL_TEXTURE_GEN_T);
	disable_blend();
	glDisable(light);
	if (depth_clamp_enabled) {glEnable(GL_DEPTH_CLAMP);}
}


void compute_brightness() {

	brightness = 0.8 + 0.2*light_factor;
	if (!have_sun) brightness *= 0.25;
	if (is_cloudy) brightness *= 0.5;
	float const sun_bright(0.5 + 0.5*max(0.0f, sun_pos.z/sun_pos.mag()));
	float const moon_bright(combined_gu ? 0.1 : 0.3*(0.5 + 0.5*max(0.0f, moon_pos.z/moon_pos.mag())));
	
	if (light_factor >= 0.6) {
		brightness *= sun_bright;
	}
	else if (light_factor <= 0.4) {
		brightness *= moon_bright;
	}
	else {
		brightness *= 5.0*((light_factor - 0.4)*sun_bright + (0.6 - light_factor)*moon_bright);
	}
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
		set_color(color2);
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
		if (status && sphere_in_camera_view(pos, radius, 0)) {
			draw_part(pos, radius, color, qbd);
		}
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
		
		if (!outside_scene && last_cid >= 0) {
			assert(last_cid < (int)coll_objects.size());
			known_coll = coll_objects[last_cid].line_intersect(p, lpos);
		}
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
	qbd.add_billboard(p, camera, up_vector, c, 4.0*r, 4.0*r, tex_range_t(), MIN_PARTICLE_FILL);
}


colorRGBA fire::get_fire_color() const { // unused

	float const alpha(rand_uniform(max(0.3, (0.9 + 0.1*heat)), min(0.9, (0.8 + 0.2*heat))));
	return colorRGBA(1.0, 0.4*heat, max(0.0f, 1.2f*(heat-1.0f)), alpha);
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


void decal_obj::draw(quad_batch_draw &qbd) const {

	assert(status);
	point const cur_pos(get_pos());
	if (dot_product_ptv(orient, cur_pos, get_camera_pos()) > 0.0) return; // back face culling
	float const alpha_val(get_alpha());
	if (!dist_less_than(cur_pos, get_camera_pos(), max(window_width, window_height)*radius*alpha_val)) return; // distance culling
	colorRGBA draw_color(color);
	draw_color.alpha = alpha_val;
	vector3d upv(orient.y, orient.z, orient.x); // swap the xyz values to get an orthogonal vector
	if (rot_angle != 0.0) {rotate_vector3d(orient, rot_angle, upv);}
	// move slightly away from the object to blend properly with cracks
	qbd.add_billboard((cur_pos + DECAL_OFFSET*orient), (cur_pos + orient), upv, draw_color, radius, radius, tex_range);
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

	if (bubbles.empty()) return;
	glEnable(GL_CULL_FACE);
	enable_blend();
	set_color(WATER_C);
	bool const set_liquid_color(world_mode == WMODE_GROUND);
	draw_objects(bubbles, set_liquid_color);
	bubble_pld.draw_and_clear();
	disable_blend();
	glDisable(GL_CULL_FACE);
}


void draw_part_clouds(vector<particle_cloud> const &pc, colorRGBA const &color, bool zoomed) {

	enable_flares(color, zoomed); // color will be set per object
	//select_multitex(CLOUD_TEX, 1);
	glAlphaFunc(GL_GREATER, 0.01);
	glEnable(GL_ALPHA_TEST); // makes it faster
	quad_batch_draw qbd;
	draw_objects(pc, qbd);
	qbd.draw();
	glDisable(GL_ALPHA_TEST);
	disable_flares();
	//set_active_texture(0);
}


void water_particle_manager::draw() const {

	if (parts.empty()) return;
	// use point sprites?
	// calculate normal in a shader?
	point const camera(get_camera_pos());
	vector<vert_norm_color> verts;
	verts.resize(parts.size());

	for (unsigned i = 0; i < parts.size(); ++i) {
		vector3d const p2c(camera - parts[i].p);
		verts[i] = vert_norm_color(parts[i].p, p2c.get_norm(), parts[i].c.c); // normal faces camera
		verts[i].c[3] *= min(1.0, 2.0/p2c.mag());
	}
	glEnable(GL_COLOR_MATERIAL);
	glPointSize(2.0);
	enable_blend();
	draw_verts(verts, GL_POINTS);
	disable_blend();
	glPointSize(1.0);
	glDisable(GL_COLOR_MATERIAL);
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


void create_and_draw_cracks() { // adds to beams

	if (decals.empty()) return;
	vector<crack_point> cpts;  // static?
	vector<ray2d> crack_lines; // static?
	int last_cobj(-1);
	bool skip_cobj(0);
	point const camera(get_camera_pos());

	for (vector<decal_obj>::const_iterator i = decals.begin(); i != decals.end(); ++i) {
		if (i->status == 0 || !i->is_glass || i->cid < 0) continue;
		if (i->cid == last_cobj && skip_cobj)             continue;
		point const pos(i->get_pos());
		if (!dist_less_than(camera, pos, 1000*i->radius)) continue; // too far away
		assert((unsigned)i->cid < coll_objects.size());
		coll_obj const &cobj(coll_objects[i->cid]);
		skip_cobj = (cobj.status != COLL_STATIC || cobj.type != COLL_CUBE || !camera_pdu.cube_visible(cobj) || cobj.is_occluded_from_camera());
		last_cobj = i->cid;
		if (skip_cobj) continue;
		int const face(cobj.closest_face(pos)), dim(face >> 1), dir(face & 1);
		vector3d dpos(all_zeros);
		
		if ((pos[dim] - camera[dim] < 0) ^ dir) { // back facing - render the crack on the other side of the glass
			dpos = 2*(cobj.get_center_pt() - i->pos);
			UNROLL_3X(dpos[i_] *= fabs(i->orient[i_]);)
		}
		cpts.push_back(crack_point(pos+dpos, i->pos+dpos, i->cid, face, i->time, i->get_alpha(), i->color));
	}
	stable_sort(cpts.begin(), cpts.end());

	for (unsigned i = 0; i < cpts.size();) {
		unsigned const s(i);
		for (++i; i < cpts.size() && cpts[i].cid == cpts[s].cid && cpts[i].face == cpts[s].face; ++i) {}
		// all cpts in [s,i) have the same {cid, face}
		crack_lines.resize(0);
		cube_t const &cube(coll_objects[cpts[s].cid]);
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

					if (attempt == 0 || dist_sq < min_dist_sq) {
						epos = p2;
						min_dist_sq = dist_sq;
					}
				} // for attempt
				beams.push_back(beam3d(0, NO_SOURCE, cpt1.pos, epos, cpt1.color, 0.05*cpt1.alpha));
				epts[n] = epos;
			} // for n
			for (unsigned n = 0; n < ncracks; ++n) {
				crack_lines.push_back(ray2d(x1, y1, epts[n][d1], epts[n][d2]));
			}
		} // for j
	} // for i
}


void draw_cracks_and_decals() {

	if (decals.empty()) return;
	create_and_draw_cracks(); // adds to beams
	map<int, quad_batch_draw> batches; // maps from {tid, is_black} to quad batches
	vector<pair<int, unsigned> > sorted_decals;

	for (obj_vector_t<decal_obj>::const_iterator i = decals.begin(); i != decals.end(); ++i) {
		if (i->status && sphere_in_camera_view(i->get_pos(), i->radius, 0)) {
			sorted_decals.push_back(make_pair(-i->time, (i - decals.begin()))); // negate time, so largest time is first
		}
	}
	if (sorted_decals.empty()) return;
	sort(sorted_decals.begin(), sorted_decals.end()); // sort by time, so that spraypaint works (later paint is drawn after/over earlier paint)

	for (unsigned i = 0; i < sorted_decals.size(); ++i) {
		decal_obj const &d(decals[sorted_decals[i].second]);
		d.draw(batches[(d.tid << 1) + (d.color == BLACK)]);
	}
	set_color(BLACK);
	glDepthMask(GL_FALSE);
	enable_blend();
	shader_t black_shader, lighting_shader, bullet_shader;

	for (map<int, quad_batch_draw>::const_iterator i = batches.begin(); i != batches.end(); ++i) {
		int const tid(i->first >> 1);
		bool const is_black(i->first & 1);

		if (tid == BULLET_D_TEX) {
			if (!bullet_shader.is_setup()) {
				// see http://cowboyprogramming.com/2007/01/05/parallax-mapped-bullet-holes/
				bullet_shader.set_prefix("#define TEXTURE_ALPHA_MASK",  1); // FS
				bullet_shader.set_prefix("#define ENABLE_PARALLAX_MAP", 1); // FS
				setup_smoke_shaders(bullet_shader, 0.05, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1); // bump maps enabled
				bullet_shader.add_uniform_float("ambient_scale", 0.0); // shader is unique, so we don't have to reset this
				bullet_shader.add_uniform_float("bump_tb_scale", -1.0); // invert the coordinate system (FIXME: something backwards?)
				bullet_shader.add_uniform_float("hole_depth", 0.2);
				bullet_shader.add_uniform_int("depth_map", 9);
				set_active_texture(5);
				select_texture(BULLET_N_TEX, 0);
				set_active_texture(9);
				select_texture(BULLET_D_TEX, 0);
				set_active_texture(0);
			}
			bullet_shader.enable();
		}
		else if (is_black) {
			if (!black_shader.is_setup()) {setup_smoke_shaders(black_shader, 0.01, 0, 1, 0, 0, 0, 1);} // no lighting
			black_shader.enable();
		}
		else {
			if (!lighting_shader.is_setup()) {
				setup_smoke_shaders(lighting_shader, 0.01, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1);
				lighting_shader.add_uniform_float("ambient_scale", 0.0);
			}
			lighting_shader.enable();
		}
		select_texture(tid, 0, 1);
		i->second.draw();
	} // for i
	disable_blend();
	glDepthMask(GL_TRUE);
	if (bullet_shader.is_setup()  ) {bullet_shader.enable();   bullet_shader.add_uniform_float  ("bump_tb_scale", 1.0);} // reset
	if (lighting_shader.is_setup()) {lighting_shader.enable(); lighting_shader.add_uniform_float("ambient_scale", 1.0);} // reset
	black_shader.end_shader();
	lighting_shader.end_shader();
	bullet_shader.end_shader();
}


void draw_smoke_and_fires() {

	if (part_clouds.empty() && fires.empty()) return; // nothing to draw
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 1, 0, 0, 0, 1);
	set_color(BLACK);

	if (!part_clouds.empty()) { // Note: just because part_clouds is nonempty doesn't mean there is any enabled smoke
		draw_part_clouds(part_clouds, WHITE, 0); // smoke: slow when a lot of smoke is up close
	}
	order_vect_t fire_order;
	get_draw_order(fires, fire_order);
	
	if (!fire_order.empty()) {
		enable_blend();
		quad_batch_draw qbd;
		select_texture(FIRE_TEX, 0);
		int last_in_smoke(-1);
		for (unsigned j = 0; j < fire_order.size(); ++j) {fires[fire_order[j].second].draw(qbd, last_in_smoke);}
		qbd.draw();
		set_std_blend_mode();
		disable_blend();
	}
	s.end_shader();
}


void add_camera_filter(colorRGBA const &color, unsigned time, int tid, unsigned ix, bool fades) {
	
	assert(ix < MAX_CFILTERS);
	if (color.alpha == 0.0) return;
	if (cfilters.size() <= ix) cfilters.resize(ix+1);
	cfilters[ix] = camera_filter(color, time, tid, fades);
}


void camera_filter::draw() {

	if (tid >= 0) select_texture(tid);
	float const zval(-1.1*perspective_nclip), tan_val(tan(perspective_fovy/TO_DEG));
	float const y(0.5*zval*tan_val), x((y*window_width)/window_height);
	colorRGBA cur_color(color);
	if (fades) {cur_color.alpha *= float(time)/float(init_time);}
	cur_color.do_glColor();
	draw_tquad(x, y, zval);
	if (tid >= 0) glDisable(GL_TEXTURE_2D);
}


void draw_camera_filters(vector<camera_filter> &cfs) {

	if (cfs.empty()) return;
	GLboolean lighting(glIsEnabled(GL_LIGHTING));
	if (lighting) glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	enable_blend();

	for (int i = (int)cfs.size()-1; i >= 0; --i) { // apply backwards
		if (cfs[i].time == 0) continue;
		cfs[i].draw();
		if ((int)cfs[i].time <= iticks) cfs[i].time = 0; else cfs[i].time -= iticks;
	}
	disable_blend();
	glEnable(GL_DEPTH_TEST);
	if (lighting) glEnable(GL_LIGHTING);
}


float const spark_t::radius = 0.0;


void spark_t::draw(quad_batch_draw &qbd) const {

	point const camera(get_camera_pos());
	qbd.add_billboard((pos + (camera - pos).get_norm()*0.02), camera, up_vector, c, 0.8*s, 0.8*s);
}


void draw_sparks() { // FIXME SHADERS: uses fixed function pipeline

	if (sparks.empty()) return;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDisable(GL_LIGHTING);
	enable_blend();
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.01);
	set_additive_blend_mode();
	select_texture(FLARE2_TEX);
	quad_batch_draw qbd;
	draw_objects(sparks, qbd);
	qbd.draw();
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glDisable(GL_ALPHA_TEST);
	set_std_blend_mode();
	disable_blend();
	set_fill_mode();
	sparks.clear();
}


void draw_projectile_effects() {

	update_blasts(); // not really an update, but needed for draw_blasts
	draw_blasts();
	draw_beams();
	draw_sparks();
	water_part_man.draw(); // not really a projectile effect, but it's drawn with them
}


void draw_env_other() {

	if (!enable_fsource) return;
	set_color(BLACK);
	draw_subdiv_sphere(flow_source, 0.05, N_SPHERE_DIV, 0, 0);
}


void draw_splash(float x, float y, float z, float size, colorRGBA color) {

	assert(size >= 0.0);
	if (DISABLE_WATER || !(display_mode & 0x04)) return;
	if (size == 0.0 || temperature <= W_FREEZE_POINT) return;
	if (size > 0.1) size = sqrt(10.0*size)/10.0;
	unsigned const num_rings(min(10U, (unsigned)ceil(size)));
	size = min(size, 0.025f);
	float radius(size);
	float const dr(0.5*size);
	point const pos(x, y, z+SMALL_NUMBER);
	unsigned const ndiv(max(3, min(N_CYL_SIDES, int(1000.0*size/max(TOLERANCE, distance_to_camera(pos))))));
	select_liquid_color(color, get_xpos(x), get_ypos(y));
	set_color(color);
	set_fill_mode();
	glPushMatrix();
	translate_to(pos);

	for (unsigned i = 0; i < num_rings; ++i) {
		draw_circle_normal((radius - 0.5*dr), radius, ndiv, 0);
		radius += dr;
	}
	glPopMatrix();
}


void draw_text(float x, float y, float z, char const *text, float tsize, bool bitmap_font) {

	//bitmap_font |= ((display_mode & 0x80) != 0);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	if (bitmap_font) {
		glRasterPos3f(x, y, z);
	}
	else {
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glPushMatrix();
		glTranslatef(x, y, z);
		uniform_scale(0.000005*tsize);
	}
	unsigned line_num(0);

	while (*text) {
		if (*text == '\n') { // newline (CR/LF)
			++line_num;

			if (bitmap_font) {
				glRasterPos3f(x, y-(0.5*line_num)/window_height, z);
			}
			else {
				glPopMatrix();
				glPushMatrix();
				glTranslatef(x, y-0.001*line_num*tsize, z);
				uniform_scale(0.000005*tsize);
			}
		}
		else {
			if (bitmap_font) {
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *text); // other fonts available
			}
			else {
				glutStrokeCharacter(GLUT_STROKE_ROMAN, *text); // GLUT_STROKE_MONO_ROMAN
			}
		}
		text++;
	}
	if (!bitmap_font) {
		glPopMatrix();
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
	}
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
}


void draw_framerate(float val) {

	char text[32];
	WHITE.do_glColor();
	sprintf(text, "%3.1f", val);
	float const ar(((float)window_width)/((float)window_height));
	draw_text(-0.011*ar, -0.011, -2.0*NEAR_CLIP, text);
}


void draw_compass_and_alt() { // and temperature

	char text[64];
	float const aspect_ratio((float)window_width/(float)window_height);
	string const dirs[8] = {"N", "NW", "W", "SW", "S", "SE", "E", "NE"};
	YELLOW.do_glColor();
	sprintf(text, "Loc: (%3.2f, %3.2f, %3.2f)", (camera_origin.x+(xoff2-xoff)*DX_VAL), (camera_origin.y+(yoff2-yoff)*DY_VAL), camera_origin.z);
	draw_text(-0.005*aspect_ratio, -0.01, -0.02, text);
	float const theta(safe_acosf(cview_dir.x)*TO_DEG);
	int const octant(int(((cview_dir.y < 0) ? (360.0 - theta) : theta)/45.0 + 22.5)&7);
	sprintf(text, "%s", dirs[octant].c_str());
	draw_text(0.005*aspect_ratio, -0.01, -0.02, text);
	sprintf(text, "Temp: %iC", int(temperature));
	draw_text(0.007*aspect_ratio, -0.01, -0.02, text);
}



