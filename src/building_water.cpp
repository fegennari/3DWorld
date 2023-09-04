// 3D World - Building Basement Water Update and Drawing
// by Frank Gennari 9/1/23

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "buildings.h"
#include "openal_wrap.h"

unsigned const MAX_SPLASHES = 40; // must agree with fragment shader code


extern int player_in_basement, animate2, display_mode;
extern unsigned room_mirror_ref_tid;
extern float fticks, CAMERA_RADIUS, water_plane_z;
extern building_t const *player_building;

void set_interior_lighting(shader_t &s, bool have_indir);
void reset_interior_lighting_and_end_shader(shader_t &s);
// from postproc_effects.cpp
void bind_frame_buffer_RGB(unsigned tu_id);
void apply_player_underwater_effect(colorRGBA const &color_mod);


class building_splash_manager_t {
	struct splash_t {
		float x, y, radius, height;
		cube_t bounds;
		splash_t(float x_, float y_, float r, float h, cube_t const &b) : x(x_), y(y_), radius(r), height(h), bounds(b) {}
		bool operator<(splash_t const &s) const {return (height < s.height);} // compare by height for min_element
		vector4d get_loc_rh_as_vec4() const {return vector4d(x, y, radius, min(height, 1.0f));}
		vector4d get_bounds_as_vec4() const {return vector4d(bounds.x1(), bounds.y1(), bounds.x2(), bounds.y2());}
	};
	vector<splash_t> splashes;
	unsigned last_size=0;
	float time=0.0;
public:
	void add_splash(point const &pos, float radius, float height, cube_t const &bounds) {
		assert(splashes.size() <= MAX_SPLASHES);

		if (!splashes.empty()) {
			splash_t &prev(splashes.back());

			if (dist_xy_less_than(pos, point(prev.x, prev.y, pos.z), 0.25*radius) && prev.radius < 2.0*radius) { // merge with previous splash (optimization)
				prev.height += height*(radius*radius/(prev.radius*prev.radius)); // add height scaled by surface area to add volumes
				prev.bounds.union_with_cube(bounds);
				return;
			}
		}
		splashes.emplace_back(pos.x, pos.y, radius, height, bounds);
		if (splashes.size() > MAX_SPLASHES) {splashes.erase(min_element(splashes.begin(), splashes.end()));} // limit size to MAX_SPLASHES
	}
	void next_frame(float ref_dist) { // floor_spacing can be used
		if (splashes.empty()) return;
		time += fticks;
		if (time > 600*TICKS_PER_SECOND) {time = 0.0;} // reset after 10 minutes to avoid FP precision problems
		float const timestep(min(fticks, 4.0f)/TICKS_PER_SECOND); // clamp fticks to 100ms
		float const exp_dist(0.25*ref_dist*timestep);

		for (splash_t &s : splashes) {
			float const prev_area(s.radius*s.radius);
			s.radius += exp_dist;
			s.height *= prev_area/(s.radius*s.radius); // volume preserving
		}
		splashes.erase(std::remove_if(splashes.begin(), splashes.end(), [](splash_t const &s) {return (s.height < 0.0005);}), splashes.end());
	}
	void set_shader_uniforms(shader_t &s) {
		assert(splashes.size() <= MAX_SPLASHES);
		unsigned const iter_end(max((unsigned)splashes.size(), last_size));
		char str[32] = {};

		for (unsigned i = 0; i < iter_end; ++i) { // add splashes for this frame
			sprintf(str, "splashes[%u].loc_rh", i);
			s.add_uniform_vector4d(str, ((i < splashes.size()) ? splashes[i].get_loc_rh_as_vec4() : vector4d())); // set unused slots to all zeros
			sprintf(str, "splashes[%u].bounds", i);
			s.add_uniform_vector4d(str, ((i < splashes.size()) ? splashes[i].get_bounds_as_vec4() : vector4d())); // set unused slots to all zeros
		}
		s.add_uniform_float("time", time/TICKS_PER_SECOND);
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
// Note: player steps call this function directly; all others go through building_t::check_for_water_splash()
void register_building_water_splash(point const &pos, float size, bool alert_zombies) { // Note: pos is in camera space
	if (player_building == nullptr) return; // shouldn't happen?
	cube_t const bounds(player_building->calc_splash_bounds(pos));
	if (bounds == cube_t()) return; // shouldn't happen?
	building_splash_manager.add_splash(pos, 0.5*CAMERA_RADIUS, size, bounds);
	if (alert_zombies) {register_building_sound((pos + get_tiled_terrain_model_xlate()), 0.5*size);} // alert zombies; convert pos to building space
#pragma omp critical(gen_sound)
	gen_sound_random_var(SOUND_SPLASH2, pos, 0.3*size, 0.9);
}
bool building_t::check_for_water_splash(point const &pos_bs, float size, bool full_room_height, bool draw_splash, bool alert_zombies) const { // Note: pos in building space
	if (this != player_building)    return 0; // only splashes for the building the player is in
	if (!water_visible_to_player()) return 0; // only if water is visible
	if (!point_in_water_area(pos_bs, full_room_height)) return 0;
	vector3d const xlate(get_tiled_terrain_model_xlate());
	register_building_water_splash((pos_bs + xlate), size, alert_zombies);

	if (draw_splash) { // Note: doesn't apply to player steps
		float const radius(0.05*min(size, 1.5f)*get_window_vspace());
		point const splash_pos(pos_bs.x, pos_bs.y, (interior->water_zval + 0.1*radius)); // slightly above the water surface
		interior->room_geom->particle_manager.add_particle(splash_pos, zero_vector, WHITE, radius, PART_EFFECT_SPLASH);
	}
	return 1;
}
bool building_t::point_in_water_area(point const &p, bool full_room_height) const {
	return (has_water() && get_water_cube(full_room_height).contains_pt(p));
}

cube_t building_t::calc_splash_bounds(point const &pos) const {
	if (!interior || !point_in_water_area(pos)) return cube_t(); // error?
	index_pair_t start, end;
	get_pgbr_wall_ix_for_pos(pos, start, end);
	vect_cube_t walls[2];

	for (unsigned d = 0; d < 2; ++d) { // copy sub-ranges of walls, since clip_ray_to_walls() doesn't take iterators; maybe it should?
		vect_cube_t const &pbgr_walls(interior->room_geom->pgbr_walls[d]);
		walls[d].insert(walls[d].end(), pbgr_walls.begin()+start.ix[d], pbgr_walls.begin()+end.ix[d]);
	}
	unsigned const NUM_RAYS = 90;
	cube_t const &extb(interior->basement_ext_bcube);
	float const ray_len(extb.dx()*extb.dx() + extb.dy()*extb.dy()); // max room diagonal
	cube_t bounds(pos, pos);

	for (unsigned n = 0; n < NUM_RAYS; ++n) {
		float const angle(TWO_PI*n/NUM_RAYS);
		point p2(pos + point(ray_len*sin(angle), ray_len*cos(angle), 0.0));
		float tmin(0.0), tmax(1.0);
		get_line_clip_xy(pos, p2, extb.d, tmin, tmax);
		p2 = pos + (p2 - pos)*tmax;
		clip_ray_to_walls(pos, p2, walls);
		bounds.union_with_pt(p2);
	} // for n
	return bounds; // zvals should be unused
}

void building_t::draw_water(vector3d const &xlate) const {
	if (!(display_mode & 0x04)) return; // water disabled
	
	if (!water_visible_to_player()) {
		if (player_in_basement < 3) {clear_building_water_splashes();} // clear if player has exited the extended basement
		return;
	}
	float const floor_spacing(get_window_vspace());
	if (animate2) {building_splash_manager.next_frame(floor_spacing);} // maybe should do this somewhere else? or update even if water isn't visible?
	point const camera_pos(get_camera_pos());

	if (camera_pos.z < interior->water_zval) { // player under the water
		apply_player_underwater_effect(colorRGBA(0.3, 0.4, 0.8)); // light blue-ish
		float const orig_water_plane_z(water_plane_z);
		water_plane_z = interior->water_zval;
		draw_underwater_particles(interior->basement_ext_bcube.z1());
		water_plane_z = orig_water_plane_z;
		return;
	}
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
	if (use_dlights) {setup_dlight_textures(s);} // must be before set_city_lighting_shader_opts()
	set_city_lighting_shader_opts(s, lights_bcube, use_dlights, use_smap, pcf_scale);
	set_interior_lighting(s, have_indir);
	float const water_depth(interior->water_zval - (interior->basement_ext_bcube.z1() + get_fc_thickness()));
	s.add_uniform_vector3d("camera_pos",  camera_pos);
	s.add_uniform_float("water_depth",    water_depth);
	s.add_uniform_float("water_atten",    1.0/floor_spacing); // attenuates to dark blue/opaque around this distance
	s.add_uniform_color("uw_atten_max",   uw_atten_max);
	s.add_uniform_color("uw_atten_scale", uw_atten_scale);
	building_splash_manager.set_shader_uniforms(s);
	bind_frame_buffer_RGB(1); // tu_id=1
	s.add_uniform_int("frame_buffer", 1);
	// Note: this must be *after* bind_frame_buffer_RGB() because it changes the texture
	if (room_mirror_ref_tid > 0) {bind_2d_texture(room_mirror_ref_tid);} else {select_texture(WHITE_TEX);}
	s.add_uniform_int("reflection_tex", 0);
	enable_blend(); // no longer needed?
	cube_t const water(get_water_cube());
	float const x1(water.x1()), y1(water.y1()), x2(water.x2()), y2(water.y2()), z(water.z2()), tx(1.0), ty(1.0);
	vector3d const &n(plus_z);
	vert_norm_tc const verts[4] = {vert_norm_tc(point(x1, y1, z), n, 0.0, 0.0), vert_norm_tc(point(x2, y1, z), n, tx, 0.0),
		                           vert_norm_tc(point(x2, y2, z), n, tx,  ty ), vert_norm_tc(point(x1, y2, z), n, 0.0, ty)};
	draw_quad_verts_as_tris(verts, 4);
	disable_blend();
	reset_interior_lighting_and_end_shader(s);
}



