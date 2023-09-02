// 3D World - Building Basement Water Update and Drawing
// by Frank Gennari 9/1/23

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "buildings.h"
#include "openal_wrap.h"

unsigned const MAX_SPLASHES = 16; // must agree with fragment shader code


extern int player_in_basement, display_mode;
extern unsigned room_mirror_ref_tid;
extern float fticks, CAMERA_RADIUS;
extern building_t const *player_building;

void set_interior_lighting(shader_t &s, bool have_indir);
void reset_interior_lighting_and_end_shader(shader_t &s);


class building_splash_manager_t {
	struct splash_t {
		float x, y, radius, height;
		splash_t(float x_, float y_, float r, float h) : x(x_), y(y_), radius(r), height(h) {}
		bool operator<(splash_t const &s) const {return (height < s.height);} // compare by height for min_element
		vector4d as_vec4() const {return vector4d(x, y, radius, height);}
	};
	vector<splash_t> splashes;
	unsigned last_size=0;
public:
	void add_splash(point const &pos, float radius, float height) {
		assert(splashes.size() <= MAX_SPLASHES);
		splashes.emplace_back(pos.x, pos.y, radius, height);
		if (splashes.size() > MAX_SPLASHES) {splashes.erase(min_element(splashes.begin(), splashes.end()));} // limit size to MAX_SPLASHES
	}
	void next_frame() {
		for (splash_t &s : splashes) {
			// TODO: update
		}
		splashes.erase(std::remove_if(splashes.begin(), splashes.end(), [](splash_t const &s) {return (s.height < 0.01);}), splashes.end());
	}
	void set_shader_uniforms(shader_t &s) {
		assert(splashes.size() <= MAX_SPLASHES);
		unsigned const iter_end(max((unsigned)splashes.size(), last_size));
		char str[24] = {};

		for (unsigned i = 0; i < iter_end; ++i) { // add splashes for this frame
			sprintf(str, "splashes[%u]", i);
			s.add_uniform_vector4d(str, ((i < splashes.size()) ? splashes[i].as_vec4() : vector4d())); // set unused slots to all zeros
		}
		last_size = splashes.size();
	}
	void clear() {splashes.clear();} // Note: last_size is not reset
};

// maybe this should be a building interior member?
// but then we would need some hook to clear this when the player exits one building so that we don't have old data when the player enters a different building,
// and we don't have the shader in register_player_not_in_building() to clear the uniforms;
// I suppose it's more efficient to only have one of these, and it's easier for recompile to not have this in buildings.h
building_splash_manager_t building_splash_manager;

void clear_building_water_splashes() {
	building_splash_manager.clear();
}
void register_building_water_splash(point const &pos, float size, bool play_sound) { // Note: pos is in camera space
	building_splash_manager.add_splash(pos, 0.5*CAMERA_RADIUS, size);
	if (!play_sound) return;
#pragma omp critical(gen_sound)
	gen_sound_random_var(SOUND_SPLASH2, pos, 0.3*size, 0.9);
}
bool building_t::check_for_water_splash(point const &pos_bs, float size, bool full_room_height, bool play_sound) const { // Note: pos in building space
	if (this != player_building)    return 0; // only splashes for the building the player is in
	if (!water_visible_to_player()) return 0; // only if water is visible
	if (!point_in_water_area(pos_bs, full_room_height)) return 0;
	vector3d const xlate(get_tiled_terrain_model_xlate());
	register_building_water_splash((pos_bs + xlate), size, play_sound);
	if (play_sound) {register_building_sound(pos_bs, 0.5*size);} // alert zombies
	return 1;
}
bool building_t::point_in_water_area(point const &p, bool full_room_height) const {
	return (has_water() && get_water_cube(full_room_height).contains_pt(p));
}

void building_t::draw_water(vector3d const &xlate) const {
	if (!(display_mode & 0x04)) return; // water disabled
	
	if (!water_visible_to_player()) {
		if (player_in_basement < 3) {clear_building_water_splashes();} // clear if player has exited the extended basement
		return;
	}
	float const floor_spacing(get_window_vspace());
	shader_t s;
	cube_t const lights_bcube(get_building_lights_bcube());
	bool const use_dlights(!lights_bcube.is_all_zeros()), have_indir(0), use_smap(1);
	float const pcf_scale = 0.2;
	s.set_prefix("#define LINEAR_DLIGHT_ATTEN", 1); // FS; improves room lighting (better light distribution vs. framerate trade-off)
	if (set_dlights_booleans(s, use_dlights, 1, 0)) {s.set_prefix("#define NO_DL_SPECULAR", 1);} // FS
	s.set_prefix(make_shader_bool_prefix("use_shadow_map", 0), 1); // FS; not using sun/moon light, so no shadow maps
	s.set_vert_shader("building_water");
	s.set_frag_shader("ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+building_water");
	s.begin_shader();
	s.add_uniform_int("reflection_tex", 0);
	if (use_dlights) {setup_dlight_textures(s);} // must be before set_city_lighting_shader_opts()
	set_city_lighting_shader_opts(s, lights_bcube, use_dlights, use_smap, pcf_scale);
	set_interior_lighting(s, have_indir);
	if (room_mirror_ref_tid > 0) {bind_2d_texture(room_mirror_ref_tid);} else {select_texture(WHITE_TEX);}
	float const water_depth(interior->water_zval - (interior->basement_ext_bcube.z1() + get_fc_thickness()));
	s.add_uniform_vector3d("camera_pos",  get_camera_pos());
	s.add_uniform_float("water_depth",    water_depth);
	s.add_uniform_float("water_atten",    1.0/floor_spacing); // attenuates to dark blue/opaque around this distance
	s.add_uniform_color("uw_atten_max",   uw_atten_max);
	s.add_uniform_color("uw_atten_scale", uw_atten_scale);
	building_splash_manager.set_shader_uniforms(s);
	enable_blend();
	cube_t const water(get_water_cube());
	float const x1(water.x1()), y1(water.y1()), x2(water.x2()), y2(water.y2()), z(water.z2()), tx(1.0), ty(1.0);
	vector3d const &n(plus_z);
	vert_norm_tc const verts[4] = {vert_norm_tc(point(x1, y1, z), n, 0.0, 0.0), vert_norm_tc(point(x2, y1, z), n, tx, 0.0),
		                           vert_norm_tc(point(x2, y2, z), n, tx,  ty ), vert_norm_tc(point(x1, y2, z), n, 0.0, ty)};
	draw_quad_verts_as_tris(verts, 4);
	disable_blend();
	reset_interior_lighting_and_end_shader(s);
}



