// 3D World - Building Basement Water Update and Drawing
// by Frank Gennari 9/1/23

#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"
#include "buildings.h"
#include "openal_wrap.h"

unsigned const MAX_SPLASHES = 40; // must agree with fragment shader code
float const BUBBLE_VELOCITY = 0.0002;


extern int player_in_basement, player_in_water, animate2, display_mode, color_buffer_frame;
extern unsigned room_mirror_ref_tid;
extern float fticks, CAMERA_RADIUS, water_plane_z;
extern building_t const *player_building;

float get_player_oxygen();
void setup_building_draw_shader_post(shader_t &s, bool have_indir);
void reset_interior_lighting_and_end_shader(shader_t &s);
// from postproc_effects.cpp
void bind_frame_buffer_RGB(unsigned tu_id);
void apply_player_underwater_effect(colorRGBA const &color_mod, float intensity=1.0);
void add_postproc_underwater_fog(float atten_scale, float max_uw_dist, float mud_amt);
void bind_depth_buffer(unsigned tu_id);


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
	void next_frame(float ref_dist, bool is_pool) { // floor_spacing can be used
		if (splashes.empty()) return;
		time += fticks;
		if (time > 600*TICKS_PER_SECOND) {time = 0.0;} // reset after 10 minutes to avoid FP precision problems
		float const timestep(min(fticks, 4.0f)/TICKS_PER_SECOND); // clamp fticks to 100ms
		float const exp_dist(0.25*ref_dist*timestep);

		for (splash_t &s : splashes) {
			float const prev_area(s.radius*s.radius);
			s.radius += exp_dist;
			float height_change(prev_area/(s.radius*s.radius)); // volume preserving
			if (is_pool) {height_change = 0.25 + 0.75*height_change;} // slower falloff in pool due to reflections off sides that aren't accounted for
			s.height *= height_change;
		} // for s
		splashes.erase(std::remove_if(splashes.begin(), splashes.end(), [](splash_t const &s) {return (s.height < 0.0005);}), splashes.end());
	}
	void set_shader_uniforms(shader_t &s, bool is_pool) {
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
		s.add_uniform_float("ripple_freq", (is_pool ? 20.0 : 10.0)); // higher frequency ripples in pools
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
	point const pos_bs(pos - get_tiled_terrain_model_xlate());
	cube_t const bounds(player_building->calc_splash_bounds(pos_bs));
	if (bounds == cube_t()) return; // shouldn't happen?
	building_splash_manager.add_splash(pos_bs, 0.5*CAMERA_RADIUS, size, bounds);
	if (alert_zombies) {register_building_sound(pos_bs, 0.5*size);} // alert zombies; convert pos to building space
#pragma omp critical(gen_sound)
	gen_sound_random_var(SOUND_SPLASH2, pos, 0.3*size, 0.9);
}
bool building_t::check_for_water_splash(point const &pos_bs, float size, bool full_room_height, bool draw_splash, bool alert_zombies) const { // Note: pos in building space
	if (player_in_water == 2)       return 0; // no splash sound if the player is underwater
	if (this != player_building)    return 0; // only splashes for the building the player is in
	if (!water_visible_to_player()) return 0; // only if water is visible
	if (!point_in_water_area(pos_bs, full_room_height)) return 0;
	vector3d const xlate(get_tiled_terrain_model_xlate());
	register_building_water_splash((pos_bs + xlate), size, alert_zombies);

	if (draw_splash && has_room_geom()) { // Note: doesn't apply to player steps
		float const radius(0.1*min(size, 1.5f)*get_window_vspace());
		point const splash_pos(pos_bs.x, pos_bs.y, (interior->water_zval + 1.5*radius)); // above the water surface (final radius is 2x original)
		interior->room_geom->particle_manager.add_particle(splash_pos, zero_vector, WHITE, radius, PART_EFFECT_SPLASH);
	}
	return 1;
}
cube_t building_t::calc_splash_bounds(point const &pos) const {
	if (!interior || !point_in_water_area(pos)) return cube_t(); // error?
	if (has_pool()) return interior->pool; // no walls
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

bool building_t::point_in_water_area(point const &p, bool full_room_height) const {
	return (has_water() && get_water_cube(full_room_height).contains_pt(p));
}
bool building_t::set_float_height(point &pos, float radius, float ceil_zval, float density) const { // density in (0.0, 1.0]
	assert(density > 0.0);
	if (density >= 1.0) return 0; // sinks
	if (!point_in_water_area((pos - radius*plus_z), 0)) return 0; // test bottom point; full_room_height=0
	max_eq(pos.z, (interior->water_zval + radius*(1.0f - 2.0f*density))); // floats on the water
	if (radius > 0.0 && !has_pool()) {min_eq(pos.z, ceil_zval - radius);} // if water level is high, and this is for backrooms (not a pool), keep below the ceiling
	return 1;
}
float building_t::get_floor_below_water_level() const {
	assert(has_water());
	unsigned const floor_ix(get_ext_basement_floor_ix(interior->water_zval));
	return interior->basement_ext_bcube.z1() + floor_ix*get_window_vspace();
}
cube_t building_t::get_water_cube(bool full_room_height) const {
	if (!has_water()) return cube_t(); // no water; error?

	if (has_pool()) {
		cube_t water(interior->pool);
		water.z2() = (full_room_height ? (water.z2() + get_floor_ceil_gap()) : interior->water_zval);
		return water;
	}
	assert(has_ext_basement()); // backrooms
	cube_t water(interior->basement_ext_bcube);
	if (full_room_height) {water.z2() = get_floor_below_water_level() + get_window_vspace();} // floor above water level
	else {water.z2() = interior->water_zval;}
	return water;
}
bool building_t::water_visible_to_player() const {
	if (!has_water()) return 0;
	vector3d const xlate(get_tiled_terrain_model_xlate());
	point const camera_bs(camera_pdu.pos - xlate);
	if (point_in_water_area(camera_bs)) return 1; // definitely visible
	if (!point_in_extended_basement_not_basement(camera_bs)) return 0;
	float const floor_spacing(get_window_vspace()), floor_below(get_floor_below_water_level()), floor_above(floor_below + floor_spacing);
	if (camera_bs.z > floor_above + floor_spacing)     return 0; // player not on the floor with water or the floor above (in case water is visible through stairs)
	if (!is_rot_cube_visible(get_water_cube(), xlate)) return 0;

	if (interior->has_backrooms) { // backrooms water
		for (stairwell_t const &s : interior->stairwells) { // check stairs visibility; stairs should be straight
			if (s.z1() > interior->water_zval) continue; // above the water level
			if (s.z2() < floor_above)          continue; // stairs don't go up to the floor the player is on
			if (!interior->basement_ext_bcube.contains_cube(s))          continue; // not extended basement stairs
			if (!s.closest_dist_less_than(camera_bs, 5.0*floor_spacing)) continue; // too far away
			cube_t floor_cut(s);
			set_cube_zvals(floor_cut, floor_above, floor_above+get_fc_thickness());
			if (is_rot_cube_visible(floor_cut, xlate)) return 1;
		} // for s
	}
	if (has_pool()) { // pool water
		room_t const &room(get_room(interior->pool.room_ix));
		if (room.contains_pt_exp(camera_bs, 2.0*get_wall_thickness())) return 1; // player in room with the pool, or in its doorway
		// check if pool is visible through a doorway
		cube_t pool_surface(interior->pool);
		pool_surface.z1() = interior->water_zval;

		for (door_stack_t const &ds : interior->door_stacks) {
			if (!ds.is_connected_to_room(interior->pool.room_ix)) continue;
			assert(ds.num_doors == 1); // must be a single door
			door_t const &door(get_door(ds.first_door_ix));
			if (door.open_amt == 0) continue; // fully closed
			if (!is_rot_cube_visible(door.get_true_bcube(), xlate)) continue;
			if (is_cube_visible_through_door(camera_bs, pool_surface, door)) return 1;
		} // for ds
	}
	return 0;
}

void building_t::draw_water(vector3d const &xlate) const {
	if (!(display_mode & 0x04)) return; // water disabled
	
	if (!water_visible_to_player()) {
		if (player_in_basement < 3) {clear_building_water_splashes();} // clear if player has exited the extended basement
		return;
	}
	cube_t const water(get_water_cube());
	bool const is_pool(has_pool());
	float const floor_spacing(get_window_vspace()), atten_scale((is_pool ? 0.3 : 1.0)/floor_spacing), water_z1(water.z1());
	if (animate2) {building_splash_manager.next_frame(floor_spacing, is_pool);} // maybe should do this somewhere else? or update even if water isn't visible?
	point const camera_pos(get_camera_pos());
	float const mud_amt(is_pool ? 0.0 : 0.5);

	if (camera_pos.z < interior->water_zval) { // player under the water; could also check (player_in_water == 2)
		point const camera_bs(camera_pos - get_tiled_terrain_model_xlate());

		if (animate2 && has_room_geom()) { // add bubbles
			static float next_bubble_time(0.0);

			if (tfticks > next_bubble_time) { // time for a bubble
				static rand_gen_t rgen;
				point bubble_pos(camera_bs);
				bubble_pos += 0.5*CAMERA_RADIUS*cview_dir; // in front of the player
				bubble_pos += 0.1*CAMERA_RADIUS*rgen.signed_rand_vector_spherical(); // randomize a bit
				float const radius(0.05*CAMERA_RADIUS*rgen.rand_uniform(0.5, 1.0));
				interior->room_geom->particle_manager.add_particle(bubble_pos, BUBBLE_VELOCITY*plus_z, WHITE, radius, PART_EFFECT_BUBBLE);
				if (rgen.rand_bool()) {gen_sound_thread_safe(SOUND_BUBBLE, camera_pos, 0.25);}
				next_bubble_time = tfticks + rgen.rand_uniform(0.2, 0.4)*TICKS_PER_SECOND; // random spacing in time
			}
		}
		float const oxygen(get_player_oxygen()), intensity(1.0 + 20.0*max(0.0f, (0.2f - oxygen))); // intensity increases when low on oxygen
		// compute max dist, for approximate use in view ray clipping to the water surface
		float max_uw_dist(0.0);
		point pts[8];
		unsigned const npts(get_cube_corners(water.d, pts)); // get all corners; we could use visible corners, but then there would be a pop when a corner becomes visible
		for (unsigned n = 0; n < npts; ++n) {max_eq(max_uw_dist, p2p_dist(camera_bs, pts[n]));}
		colorRGBA const pool_color(0.7, 0.8, 1.0), clear_color(0.4, 0.6, 1.0), mud_color(1.0, 0.6, 0.33);
		colorRGBA uw_color(is_pool ? pool_color : blend_color(mud_color, clear_color, mud_amt, 0));
		apply_player_underwater_effect(uw_color*min(1.0, 10.0*oxygen), intensity); // fade to black when oxygen is low
		add_postproc_underwater_fog((is_pool ? 2.0 : 1.0)*WATER_COL_ATTEN*atten_scale, max_uw_dist, mud_amt);
		bool const is_lit(is_room_lit(get_room_containing_pt(camera_bs), camera_bs.z));
		colorRGBA const base_color(is_lit ? WHITE : DK_GRAY);
		float const orig_water_plane_z(water_plane_z);
		water_plane_z = interior->water_zval;
		draw_underwater_particles(water_z1, base_color);
		water_plane_z = orig_water_plane_z;
		return;
	}
	shader_t s;
	float water_depth(0.0);

	if (is_pool) {water_depth = interior->water_zval - interior->pool.z1();}
	else { // backrooms water
		unsigned const camera_floor(unsigned((min(camera_pos.z, interior->water_zval) - water_z1)/floor_spacing)); // handle player on floor above water
		water_depth = interior->water_zval - (water_z1 + get_fc_thickness() + camera_floor*floor_spacing); // for the player's floor
		min_eq(water_depth, get_floor_ceil_gap()); // lights are on every floor, so optical depth can't be more than the distance between the floor and the lights above it
	}
	cube_t const lights_bcube(get_building_lights_bcube());
	bool const use_dlights(!lights_bcube.is_all_zeros()), have_indir(0), use_smap(1); // indir lighting has little effect and is difficult to setup
	float const pcf_scale = 0.2;
	s.set_prefix("#define LINEAR_DLIGHT_ATTEN", 1); // FS; improves room lighting (better light distribution vs. framerate trade-off)
	if (set_dlights_booleans(s, use_dlights, 1, 0)) {s.set_prefix("#define NO_DL_SPECULAR", 1);} // FS
	s.set_prefix(make_shader_bool_prefix("use_shadow_map", 0), 1); // FS; not using sun/moon light, so no shadow maps
	s.set_vert_shader("building_water");
	s.set_frag_shader("ads_lighting.part*+shadow_map.part*+dynamic_lighting.part*+depth_utils.part+building_water");
	s.begin_shader();
	if (use_dlights) {setup_dlight_textures(s);} // must be before set_city_lighting_shader_opts()
	set_city_lighting_shader_opts(s, lights_bcube, use_dlights, use_smap, pcf_scale);
	setup_building_draw_shader_post(s, have_indir);
	s.add_uniform_vector3d("camera_pos", camera_pos);
	s.add_uniform_float("water_depth",   water_depth);
	s.add_uniform_float("foam_scale",    min(1.0f, 0.1f*floor_spacing/water_depth)); // higher with shallow water, lower with deep water
	setup_shader_underwater_atten(s, atten_scale, mud_amt); // attenuates to dark blue (or brown for mud)/opaque around this distance
	building_splash_manager.set_shader_uniforms(s, is_pool);
	bind_frame_buffer_RGB(8); // tu_id=8
	s.add_uniform_int("frame_buffer", 8); // tu_id=8
	bind_depth_buffer( 9); // tu_id=9
	setup_depth_tex(s, 9); // tu_id=9
	// Note: this must be *after* bind_frame_buffer_RGB() and bind_depth_buffer() because it changes the texture
	if (room_mirror_ref_tid > 0) {bind_texture_tu(room_mirror_ref_tid, 0);} else {select_texture(WHITE_TEX);}
	s.add_uniform_int("reflection_tex", 0);
	enable_blend(); // no longer needed?
	float const x1(water.x1()), y1(water.y1()), x2(water.x2()), y2(water.y2()), z(water.z2()), tx(1.0), ty(1.0);
	vector3d const &n(plus_z);
	vert_norm_tc const verts[4] = {vert_norm_tc(point(x1, y1, z), n, 0.0, 0.0), vert_norm_tc(point(x2, y1, z), n, tx, 0.0),
		                           vert_norm_tc(point(x2, y2, z), n, tx,  ty ), vert_norm_tc(point(x1, y2, z), n, 0.0, ty)};
	draw_quad_verts_as_tris(verts, 4);
	disable_blend();
	reset_interior_lighting_and_end_shader(s);
	color_buffer_frame = 0; // reset to invalidate buffer
}



