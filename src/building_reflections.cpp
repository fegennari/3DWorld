// 3D World - Building Mirror Reflections
// by Frank Gennari 9/29/2020

#include "function_registry.h"
#include "buildings.h"
#include "shaders.h"

bool disable_city_shadow_maps(0), mirror_in_ext_basement(0);
unsigned room_mirror_ref_tid(0);
room_object_t cur_room_mirror;
shader_t reflection_shader;

extern int display_mode, window_width, window_height;
extern float CAMERA_RADIUS;
extern vector4d clip_plane;
extern building_t const *player_building;

colorRGBA get_clear_color();
cube_t get_mirror_surface(room_object_t const &c);

bool is_mirror(room_object_t const &obj) {return (obj.type == TYPE_MIRROR || obj.type == TYPE_DRESS_MIR);}

bool cube_visible_in_building_mirror_reflection(cube_t const &c) {
	if (!is_mirror(cur_room_mirror)) return 0;
	pos_dir_up pdu(camera_pdu);
	// below code is duplicated with create_mirror_reflection_if_needed()/draw_scene_for_building_reflection(), but not easy to factor out and share
	bool const is_open(cur_room_mirror.is_open()), dim(cur_room_mirror.dim ^ is_open), dir(is_open ? 1 : cur_room_mirror.dir);
	float const reflect_plane(is_open ? get_mirror_surface(cur_room_mirror).d[dim][1] : cur_room_mirror.d[dim][dir]);
	float const reflect_plane_xf(reflect_plane + get_tiled_terrain_model_xlate()[dim]);
	pdu.apply_dim_mirror(dim, reflect_plane_xf);
	return pdu.cube_visible(c);
}

void draw_scene_for_building_reflection(unsigned &ref_tid, unsigned dim, bool dir, float reflect_plane,
	bool is_house, bool interior_room, bool draw_exterior, bool is_extb, bool is_water, cube_t const &mirror)
{
	int reflection_pass(REF_PASS_ENABLED);
	if ( is_house     ) {reflection_pass |= REF_PASS_HOUSE   ;} // unused
	if ( interior_room) {reflection_pass |= REF_PASS_INTERIOR;}
	if ( is_water     ) {reflection_pass |= REF_PASS_WATER   ;}
	if ( is_extb      ) {reflection_pass |= REF_PASS_EXTB    ;}
	if (!draw_exterior) {reflection_pass |= REF_PASS_INT_ONLY;}
	unsigned const txsize(window_width), tysize(window_height); // full resolution
	vector3d const xlate(get_tiled_terrain_model_xlate());
	float const reflect_plane_xf(reflect_plane + xlate[dim]), reflect_sign(dir ? -1.0 : 1.0);
	point const old_camera_pos(camera_pos);
	pos_dir_up const old_camera_pdu(camera_pdu); // reflect camera frustum used for VFC
	camera_pdu.apply_dim_mirror(dim, reflect_plane_xf); // setup reflected camera frustum
	//camera_pos = camera_pdu.pos; // this can move the camera outside the building, so we can't use it
	pos_dir_up const refl_camera_pdu(camera_pdu);
	clip_plane      = vector4d();
	clip_plane[dim] = -reflect_sign;
	clip_plane.w    =  reflect_sign*reflect_plane;
	// Note: it may be more efficient to use an FBO here, but we would need both a color attachment (room_mirror_ref_tid) and a depth attachment (and stencil buffer?)
	// Note: clearing the buffers at this point in the control flow will discard some geometry that has already been drawn such as the sky,
	//       but these generally arent't visible from within the room containing the mirror anyway
	colorRGBA const orig_clear_color(get_clear_color());
	if (is_water) {glClearColor_rgba(GRAY);} // water reflections distort the UV and can go outside the drawn texture, so use a color that blends better than light blue
	setup_viewport_and_proj_matrix(txsize, tysize); // and clear
	if (is_water) {glClearColor_rgba(orig_clear_color);} // restore clear color
	apply_dim_mirror(dim, reflect_plane_xf); // setup mirror transform
	camera_pdu = refl_camera_pdu; // reset reflected PDU
	// draw the mirror area in the stencil buffer
	setup_stencil_buffer_write();
	glStencilOpSeparate(GL_BACK,  GL_KEEP, GL_KEEP, GL_KEEP); // ignore back faces
	glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_KEEP, GL_INCR); // mark stencil on front faces
	shader_t s;
	s.begin_shadow_map_shader(); // should be okay for stencil
	draw_simple_cube((mirror + xlate), 0); // draw translated mirror
	s.end_shader();
	end_stencil_write();
	// enable stencil test for drawing building interiors as an optimization
	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_NOTEQUAL, 0, ~0U); // keep if stencil bit has been set by the mirror draw
	glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
	// draw reflected geometry
	glEnable(GL_CLIP_DISTANCE0);
	draw_buildings(0, reflection_pass, xlate);
	glDisable(GL_CLIP_DISTANCE0);

	if (draw_exterior) {
		disable_city_shadow_maps = 1; // shadows don't work due to the mirror transform and are disabled for both the terrain and the city roads/objects
		if (world_mode == WMODE_INF_TERRAIN) {draw_city_roads(1, xlate);} // opaque only
		draw_tiled_terrain(2); // reflection_pass=2
		draw_building_lights(xlate);
		//draw_tiled_terrain_clouds(1); // clouds are unlikely to be reflected in mirrors so probably don't need to be drawn
		disable_city_shadow_maps = 0;
	}
	glDisable(GL_STENCIL_TEST);
	// write reflection to a texture and reset the state
	setup_reflection_texture(ref_tid, txsize, tysize);
	render_to_texture(ref_tid, txsize, tysize); // render reflection to texture
	restore_matrices_and_clear(); // reset state
	camera_pos = old_camera_pos;
	camera_pdu = old_camera_pdu; // restore camera_pdu
	clip_plane = vector4d(); // reset to disable
	if (draw_exterior) {draw_cloud_planes(0.0, 0, 1, 1);} // redraw cloud planes since they got overwritten; terrain_zmin=0 (use prev)
}

void create_mirror_reflection_if_needed() {
	if (player_building == nullptr) return; // what about pool visible from connected building extended basement?
	
	if (player_building->water_visible_to_player()) { // draw water plane reflection
		if (get_camera_pos().z < player_building->interior->water_zval + player_building->get_window_vspace()) { // only if the player is on the same floor as the water
			cube_t water_cube(player_building->get_water_cube(0));
			water_cube.z1() = water_cube.z2(); // top surface only
			mirror_in_ext_basement = 1; // required when extended basement goes outside the building's tile
			draw_scene_for_building_reflection(room_mirror_ref_tid, 2, 1, water_cube.z2(), 0, 1, 0, 1, 1, water_cube); // +z, not house, interior, basement, no exterior
			return;
		}
	}
	if (!is_mirror(cur_room_mirror)) return; // not enabled
	bool const interior_room(cur_room_mirror.is_interior()), is_house(cur_room_mirror.is_house()), is_open(cur_room_mirror.is_open());
	// assumes mirror is not facing the doorway to a room with a window; assumes cube-shaped office buildings always use opaque glass block windows
	bool const can_see_out_windows((is_house || !player_building->is_cube()) && !interior_room && player_building->has_int_windows());
	bool const is_extb(player_building->point_in_extended_basement_not_basement(cur_room_mirror.get_cube_center()));
	bool const dim(cur_room_mirror.dim ^ is_open), dir(is_open ? 1 : cur_room_mirror.dir); // always opens in +dir
	cube_t const mirror_surface(get_mirror_surface(cur_room_mirror));
	float const reflect_plane(is_open ? mirror_surface.d[dim][1] : cur_room_mirror.d[dim][dir]);
	draw_scene_for_building_reflection(room_mirror_ref_tid, dim, dir, reflect_plane, is_house, interior_room, can_see_out_windows, is_extb, 0, mirror_surface);
	cur_room_mirror = room_object_t(); // reset for next frame
}

bool building_t::line_intersect_walls(point const &p1, point const &p2, bool same_room) const {
	cube_t const line_bcube(p1, p2);
	
	if (!same_room) { // no need to check walls if both points are in the same room
		for (unsigned d = 0; d < 2; ++d) {
			if (line_int_cubes(p1, p2, interior->walls[d], line_bcube)) return 1;
		}
	}
	if (has_room_geom() && (is_pos_in_pg_or_backrooms(p1) || is_pos_in_pg_or_backrooms(p2))) {
		for (unsigned d = 0; d < 2; ++d) {
			if (line_int_cubes(p1, p2, interior->room_geom->pgbr_walls[d], line_bcube)) return 1;
		}
	}
	return 0;
}
bool building_t::is_cube_face_visible_from_pt(cube_t const &c, point const &p, unsigned dim, bool dir, bool same_room) const { // approximate
	if (same_room && !is_pos_in_pg_or_backrooms(p)) return 1; // skip intersection tests
	assert(dim < 2); // X or Y only
	unsigned const steps(21), d1(1-dim);
	float const delta(c.get_sz_dim(d1)/(steps-1));
	point cpt;
	cpt.z    = c.zc(); // no need to test all zvals since walls span the entire room height
	cpt[dim] = c.d[dim][dir]; // face plane

	for (unsigned i = 0; i < steps; ++i) {
		cpt[d1] = c.d[d1][0] + i*delta;
		if (!line_intersect_walls(p, cpt)) return 1; // this point is visible
	}
	return 0;
}

bool building_t::find_mirror_in_room(unsigned room_id, vector3d const &xlate, float &dmin_sq, bool same_room) const { // in view of the player
	assert(has_room_geom());
	point camera_bs(camera_pdu.pos - xlate);
	maybe_inv_rotate_point(camera_bs); // rotate camera pos into building space
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	float const camera_z1(camera_bs.z - CAMERA_RADIUS), camera_z2(camera_bs.z + CAMERA_RADIUS);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) { // see if that room contains a mirror
		if (i->room_id != room_id || !is_mirror(*i))      continue; // wrong room, or not a mirror
		if (i->z1() > camera_z2   || i->z2() < camera_z1) continue; // wrong floor
		// Note: we could probably return 0 rather than continuing after this point, but that may change if rooms with multiple mirrors are enabled
		if (((camera_bs[i->dim] - i->get_center_dim(i->dim)) < 0.0f) ^ i->dir ^ 1) continue; // back facing
		if (!camera_pdu.cube_visible(*i + xlate)) continue; // VFC
		if (!is_cube_face_visible_from_pt(*i, camera_bs, i->dim, i->dir, same_room)) continue; // visibility test (slow)
		point const center(i->get_cube_center());
		float const dsq(p2p_dist_sq(camera_bs, center));
		
		if (dmin_sq == 0.0 || dsq < dmin_sq) {
			room_t const &room(get_room(room_id));

			if (room.is_backrooms()) { // backrooms bathroom mirror
				// check for closed door line intersection; if found, mirror is not visible
				assert(interior->ext_basement_door_stack_ix >= 0 && (unsigned)interior->ext_basement_door_stack_ix < interior->door_stacks.size());
				bool found_closed_door(0);

				for (auto ds = interior->door_stacks.begin()+interior->ext_basement_door_stack_ix; ds != interior->door_stacks.end(); ++ds) {
					if (ds->z1() > i->z2() || ds->z2() < i->z1()) continue; // wrong floor
					assert(ds->num_doors == 1); // must be a single door stack
					door_t const &door(get_door(ds->first_door_ix));
					if (door.open_amt > 0.0) continue; // open, skip
					if (is_cube_visible_through_door(camera_bs, *i, door)) {found_closed_door = 1; break;}
				}
				if (found_closed_door) continue;
			}
			dmin_sq = dsq;
			cur_room_mirror        = *i;
			mirror_in_ext_basement = room.is_ext_basement();
		}
	} // for i
	return (dmin_sq > 0.0);
}

bool building_t::find_mirror_needing_reflection(vector3d const &xlate) const {
	if (!has_room_geom()) return 0; // can't have mirrors; maybe interior wasn't generated yet
	if (is_rotated())     return 0; // mirrors don't yet work in rotated buildings, so disable for now
	point const camera_bs(camera_pdu.pos - xlate);
	vector<point> points;
	if (!check_point_or_cylin_contained(camera_bs, 0.0, points, 0, 1, 0)) return 0; // camera not in the building; inc_attic=0, inc_ext_basement=1, inc_roof_acc=0
	int camera_room_ix(-1);
	float dmin_sq(0.0);
	
	// find room containing the camera; note that this applies to the entire backrooms, since it's one room
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (!r->contains_pt(camera_bs)) continue; // not the room the camera is in
		unsigned const room_ix(r - interior->rooms.begin());
		camera_room_ix = room_ix;
		if (!r->has_mirror) continue; // no mirror in this room stack
		if (find_mirror_in_room((room_ix & 255), xlate, dmin_sq, 1)) return 1; // same_room=1
	} // for r
	if (camera_room_ix < 0) return 0; // camera not in a room
	if (all_room_int_doors_closed(camera_room_ix, camera_bs.z)) return 0; // camera in room with doors closed, no other mirrors are visible
	room_t const &camera_room(get_room(camera_room_ix));
	cube_t search_area(camera_room);
	search_area.expand_by_xy(2.0*get_wall_thickness()); // expand so that it overlaps adjacent rooms
	bool found(0);

	// not found, look for an adjacent room or connecting hallway and select closest mirror; what about mirrors visible from more than one room away?
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		unsigned const room_ix(r - interior->rooms.begin());
		if ((int)room_ix == camera_room_ix) continue;

		if (camera_room.is_apt_or_hotel_room() && r->unit_id == camera_room.unit_id) { // same apartment or hotel room unit
			if (all_room_int_doors_closed(room_ix, camera_bs.z)) continue; // in room with doors closed, not visible
		}
		else {
			if (!r->intersects(search_area)) continue; // wrong room

			if (camera_room.is_hallway) { // special optimization logic for hallways (generally for office buildings, but can apply to houses as well)
				cube_t r_exp(*r);
				bool const short_dim(camera_room.dy() < camera_room.dx());
				r_exp.expand_by_xy(camera_room.get_sz_dim(short_dim));
				if (!r_exp.contains_pt(camera_bs)) continue; // camera not within the hallway across from the room
			}
			if (!are_rooms_connected(*r, camera_room, camera_bs.z, 1)) continue; // no door, or door is fully closed check_open=1
		}
		found |= find_mirror_in_room((room_ix & 255), xlate, dmin_sq, 0); // same_room=0
	} // for r
	return found; // not found
}

/*static*/ bool tid_nm_pair_t::bind_reflection_shader() {
	if (room_mirror_ref_tid == 0) {select_texture(WHITE_TEX); return 0;}
	// use a custom shader that uses screen coordinates to clip the texture to the mirror bounds; inefficient (wastes texels), but simple
	bind_2d_texture(room_mirror_ref_tid);

	if (reflection_shader.is_setup()) {reflection_shader.make_current();}
	else { // setup reflection_shader
		reflection_shader.set_vert_shader("mirror_reflection");
		reflection_shader.set_frag_shader("mirror_reflection");
		reflection_shader.begin_shader();
		reflection_shader.add_uniform_int("reflection_tex", 0);
	}
	return 1;
}

