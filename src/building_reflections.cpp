// 3D World - Building Mirror Reflections
// by Frank Gennari 9/29/2020

#include "function_registry.h"
#include "buildings.h"
#include "shaders.h"

bool disable_city_shadow_maps(0);
unsigned room_mirror_ref_tid(0);
room_object_t cur_room_mirror;
shader_t reflection_shader;

extern int display_mode, window_width, window_height;
extern float CAMERA_RADIUS;
extern vector4d clip_plane;

cube_t get_mirror_surface(room_object_t const &c);

bool is_mirror(room_object_t const &obj) {return (obj.type == TYPE_MIRROR || obj.type == TYPE_DRESS_MIR);}


void draw_mirror_to_stencil_buffer(vector3d const &xlate) {
	setup_stencil_buffer_write();
	glStencilOpSeparate(GL_BACK,  GL_KEEP, GL_KEEP, GL_KEEP); // ignore back faces
	glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_KEEP, GL_INCR); // mark stencil on front faces
	shader_t s;
	s.begin_color_only_shader();
	draw_simple_cube((get_mirror_surface(cur_room_mirror) + xlate), 0); // draw translated mirror
	s.end_shader();
	end_stencil_write();
}

void create_mirror_reflection_if_needed() {
	if (!is_mirror(cur_room_mirror)) return; // not enabled
	bool const interior_room(cur_room_mirror.is_interior()), is_house(cur_room_mirror.is_house()), is_open(cur_room_mirror.is_open());
	bool const can_see_out_windows(is_house && !interior_room); // assumes mirror is not facing the doorway to a room with a window
	bool const dim(cur_room_mirror.dim ^ is_open), dir(is_open ? 1 : cur_room_mirror.dir); // always opens in +dir
	int const reflection_pass(is_house ? 3 : (interior_room ? 2 : 1));
	vector3d const xlate(get_tiled_terrain_model_xlate());
	float const reflect_plane(is_open ? get_mirror_surface(cur_room_mirror).d[dim][1] : cur_room_mirror.d[dim][dir]);
	float const reflect_plane_xf(reflect_plane + xlate[dim]), reflect_sign(dir ? -1.0 : 1.0);
	clip_plane      = vector4d();
	clip_plane[dim] = -reflect_sign;
	clip_plane.w    = reflect_sign*reflect_plane;
	unsigned const txsize(window_width), tysize(window_height); // full resolution
	point const old_camera_pos(camera_pos);
	pos_dir_up const old_camera_pdu(camera_pdu); // reflect camera frustum used for VFC
	camera_pdu.apply_dim_mirror(dim, reflect_plane_xf); // setup reflected camera frustum
	//camera_pos = camera_pdu.pos; // this can move the camera outside the building, so we can't use it
	pos_dir_up const refl_camera_pdu(camera_pdu);
	// Note: it may be more efficient to use an FBO here, but we would need both a color attachment (room_mirror_ref_tid) and a depth attachment (and stencil buffer?)
	// Note: clearing the buffers at this point in the control flow will discard some geometry that has already been drawn such as the sky,
	//       but these generally arent't visible from within the bathroom anyway
	setup_viewport_and_proj_matrix(txsize, tysize);
	apply_dim_mirror(dim, reflect_plane_xf); // setup mirror transform
	camera_pdu = refl_camera_pdu; // reset reflected PDU
	draw_mirror_to_stencil_buffer(xlate);
	// enable stencil test for drawing building interiors as an optimization
	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_NOTEQUAL, 0, ~0U); // keep if stencil bit has been set by the mirror draw
	glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
	glEnable(GL_CLIP_DISTANCE0);
	draw_buildings(0, reflection_pass, xlate); // reflection_pass=1/2/3
	glDisable(GL_CLIP_DISTANCE0);

	if (can_see_out_windows) {
		disable_city_shadow_maps = 1; // shadows don't work due to the mirror transform and are disabled for both the terrain and the city roads/objects
		if (world_mode == WMODE_INF_TERRAIN) {draw_city_roads(1, xlate);} // opaque only
		draw_tiled_terrain(2); // reflection_pass=2
		draw_building_lights(xlate);
		//draw_tiled_terrain_clouds(1); // clouds are unlikely to be reflected in bathroom mirrors so probably don't need to be drawn
		disable_city_shadow_maps = 0;
	}
	glDisable(GL_STENCIL_TEST);
	// write reflection to a texture and reset the state
	setup_reflection_texture(room_mirror_ref_tid, txsize, tysize);
	render_to_texture(room_mirror_ref_tid, txsize, tysize); // render reflection to texture
	restore_matrices_and_clear(); // reset state
	camera_pos = old_camera_pos;
	camera_pdu = old_camera_pdu; // restore camera_pdu
	clip_plane = vector4d(); // reset to disable
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

bool building_t::find_mirror_in_room(unsigned room_id, vector3d const &xlate, bool same_room) const {
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
		cur_room_mirror = *i;
		return 1;
	} // for i
	return 0;
}

bool building_t::find_mirror_needing_reflection(vector3d const &xlate) const {
	if (!has_room_geom()) return 0; // can't have mirrors; maybe interior wasn't generated yet
	if (is_rotated())     return 0; // mirrors don't yet work in rotated buildings, so disable for now
	point const camera_bs(camera_pdu.pos - xlate);
	vector<point> points;
	if (!check_point_or_cylin_contained(camera_bs, 0.0, points, 0, 1)) return 0; // camera not in the building; inc_attic=0, inc_ext_basement=1
	float const wall_thickness(get_wall_thickness());
	
	// find room containing the camera
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (!r->has_mirror) continue; // no mirror in this room stack
		if (!r->contains_pt(camera_bs)) continue;
		if (find_mirror_in_room(((r - interior->rooms.begin()) & 255), xlate, 1)) return 1; // same_room=1
	}
	// not found, look for a connecting hallway
	for (auto h = interior->rooms.begin(); h != interior->rooms.end(); ++h) {
		if (!h->is_hallway || !h->contains_pt(camera_bs)) continue;
		cube_t hallway(*h);
		hallway.expand_by_xy(2.0*wall_thickness); // expand so that it overlaps adjacent rooms
		bool const short_dim(h->dy() < h->dx());
		float const hallway_width(h->get_sz_dim(short_dim));

		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
			if (r->is_hallway) continue; // exclude other hallways (including *h)
			if (!r->intersects(hallway)) continue; // wrong room
			cube_t r_exp(*r);
			r_exp.expand_in_dim( short_dim, hallway_width);
			r_exp.expand_in_dim(!short_dim, hallway_width); // expand in the other dim to include a bit of buffer along the hallway
			if (!r_exp.contains_pt(camera_bs)) continue; // camera not within the hallway across from the room
			if (find_mirror_in_room(((r - interior->rooms.begin()) & 255), xlate, 0)) return 1; // same_room=0
		} // for r
	} // for h
	return 0; // not found
}

bool tid_nm_pair_t::bind_reflection_shader() const {
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

