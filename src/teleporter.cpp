// 3D World - Teleporter/Jump Pad Classes and Functions
// by Frank Gennari
// 5/17/18

#include "gameplay.h"
#include "player_state.h"
#include "physics_objects.h"
#include "textures.h"
#include "openal_wrap.h"
#include "explosion.h"
#include "shaders.h"
#include "draw_utils.h"


extern bool begin_motion;
extern int game_mode, num_smileys;
extern float fticks, NEAR_CLIP, FAR_CLIP;
extern double tfticks;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];
extern player_state *sstates;
extern vector<teleporter> teleporters[3]; // static, dynamic, in-hand
extern vector<jump_pad> jump_pads;


// teleporters

bool maybe_teleport_object(point &opos, float oradius, int player_id, int type, bool small_object) {

	for (vector<teleporter>::iterator i = teleporters[0].begin(); i != teleporters[0].end(); ++i) { // static teleporters
		if (i->maybe_teleport_object(opos, oradius, player_id, small_object)) return 1; // we don't support collisions with multiple teleporters at the same time
	}
	if (type == TELEPORTER) return 0; // don't teleport teleporters
	int const group(coll_id[TELEPORTER]);
	if (group < 0) return 0;
	obj_group const &objg(obj_groups[group]);
	if (!objg.is_enabled()) return 0;
	point const orig_pos(opos);

	for (unsigned i = 0; i < objg.end_id; ++i) { // check dynamic teleporter objects
		dwobject const &obj(objg.get_obj(i));
		if (obj.disabled()) continue;
		if (obj.time < 0.1*TICKS_PER_SECOND) continue; // not yet armed (too close to shooter)
		assert(i < teleporters[1].size());

		if (teleporters[1][i].maybe_teleport_object(opos, oradius, player_id, small_object)) {
			if (player_id != NO_SOURCE) {
				sstates[player_id].last_teleporter = obj.source; // credit obj.source with any fall damage kill
				if (player_id == CAMERA_ID) {gen_delayed_sound(1.0, SOUND_FALLING, all_zeros, 2.0, 1.0, 1);} // relative to listener
				else {gen_delayed_sound(1.0, SOUND_FALLING, 0.5*(opos + orig_pos), 2.0);} // halfway down for smiley
			}
			return 1;
		}
	} // for i
	return 0;
}

void setup_dynamic_teleporters() {

	if (coll_id[TELEPORTER] < 0) return;
	obj_group const &objg(obj_groups[coll_id[TELEPORTER]]);
	if (!objg.is_enabled()) {teleporters[1].clear(); return;} // Note: no texture, free_context() doesn't need to be called
	teleporters[1].resize(objg.end_id);
	for (unsigned i = 0; i < objg.end_id; ++i) {teleporters[1][i].from_obj(objg.get_obj(i));}
}

void teleport_object(point &opos, point const &src_pos, point const &dest_pos, float oradius, int player_id) {

	bool const is_player(player_id != NO_SOURCE);
	float const gain(is_player ? 1.0 : 0.1), pitch(is_player ? 0.6 : 2.0);
	gen_sound(SOUND_POWERUP, opos, gain, pitch); // different sound?
	opos += (dest_pos - src_pos); // maintain relative distance from center (could also use opos = dest)
	gen_sound(SOUND_POWERUP, opos, gain, pitch); // different sound?
	if (is_player) {player_teleported(opos, player_id);}
	add_dynamic_light(12.0*oradius, opos, LT_BLUE);
	if (player_id != CAMERA_ID) {add_blastr(opos, plus_z, 6.0*oradius, 0.0, int(0.5*TICKS_PER_SECOND), NO_SOURCE, WHITE, BLUE, ETYPE_NUCLEAR, nullptr, 1);}
}

void teleporter::from_obj(dwobject const &obj) {

	float const prev_radius(radius);
	pos     = obj.pos;
	radius  = 2.0*object_types[obj.type].radius; // make it larger so that teleport radius is larger than coll radius
	draw_radius_scale = 1.0;
	dest    = pos;
	dest.z += 400.0*CAMERA_RADIUS; // large dz
	source  = obj.source;
	enabled = obj.enabled();
	is_portal = 0;
	if (radius != prev_radius) {setup();} // on radius change
}

bool teleporter::maybe_teleport_object(point &opos, float oradius, int player_id, bool small_object) {

	if (!enabled) return 0;
	if (!dist_less_than(pos, opos, radius+oradius)) return 0; // not close enough
	if (small_object) {opos += (dest - pos); return 1;} // no effects
	teleport_object(opos, pos, dest, oradius, player_id);
	last_used_tfticks = tfticks;
	return 1;
}


void create_portal_textures() { // static teleporters only
	for (auto i = teleporters[0].begin(); i != teleporters[0].end(); ++i) {i->create_portal_texture();}
}
void draw_teleporters() {

	if (teleporters[0].empty() && teleporters[1].empty() && teleporters[2].empty()) return;
	vpc_shader_t s;
	teleporter::shader_setup(s, 4); // RGBA noise
	s.enable();
	s.add_uniform_float("noise_scale", 1.2);
	s.set_cur_color(WHITE);
	enable_blend();

	// Note: drawn additively with no depth write, doesn't need to be back to front sorted
	for (unsigned d = ((world_mode == WMODE_GROUND) ? 0 : 1); d < 3; ++d) { // static (if ground mode), dynamic, and in-hand teleporters
		for (auto i = teleporters[d].begin(); i != teleporters[d].end(); ++i) {i->draw(s, (d > 0));}
	}
	disable_blend();
	teleporters[2].clear(); // clear in-hand teleporters (there for drawing only)
}
void free_teleporter_textures() {
	for (unsigned d = 0; d < 2; ++d) { // both static and dynamic teleporters (though only static should have textures allocated)
		for (auto i = teleporters[d].begin(); i != teleporters[d].end(); ++i) {i->free_context();}
	}
}

bool teleporter::do_portal_draw() const {return (enabled && is_portal && distance_to_camera(pos) < 60.0*radius);} // transparent, and not too far away

void teleporter::create_portal_texture() {

	if (!do_portal_draw()) return;
	if (!sphere_in_camera_view(pos, get_draw_radius(), 2)) return;
	unsigned const tex_size = 512; // shouldn't need to be too high

	if (tid == 0) { // allocate tid if needed
		setup_texture(tid, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, tex_size, tex_size, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
	}
	pos_dir_up const pdu(dest, cview_dir, plus_z, 0.5*TO_RADIANS*PERSP_ANGLE, NEAR_CLIP, FAR_CLIP); // same as camera, but 1:1 AR and at teleporter destination
	create_camera_view_texture(tid, tex_size, pdu, is_indoors);
}

void teleporter::draw(vpc_shader_t &s, bool is_dynamic) { // Note: not const or static because of tid caching for transparent case

	if (!enabled) return;
	float const ACTIVATE_DELAY = 1.0; // in seconds
	float const use_scale(CLIP_TO_01(ACTIVATE_DELAY - float(tfticks - last_used_tfticks)/TICKS_PER_SECOND));
	float const draw_radius(get_draw_radius()), light_radius(8.0*draw_radius), use_light_radius(2.0*use_scale*light_radius);
	colorRGBA const c1(blend_color(YELLOW,  RED,    use_scale, 0));
	colorRGBA const c2(blend_color(WHITE,   YELLOW, use_scale, 0));
	colorRGBA const c3(blend_color(LT_BLUE, BLUE,   use_scale, 0));

	if (use_scale > 0.0 && sphere_in_camera_view(pos, use_light_radius, 2)) {
		add_dynamic_light(use_light_radius, pos, LT_BLUE, plus_z, 1.0, nullptr, 1); // static pos
	}
	if (game_mode && begin_motion && sphere_in_camera_view(pos, light_radius, 2)) {
		colorRGBA const lt_color(blend_color(blend_color(c1, c2, fabs(sin(0.05*tfticks)), 0), c3, fabs(cos(0.07*tfticks)), 0));
		add_dynamic_light(light_radius, pos, lt_color, plus_z, 1.0, nullptr, 1); // static pos
	}
	if (sphere_in_camera_view(pos, draw_radius, 2)) { // draw pos
		s.set_uniform_color(s.c1i_loc, c1);
		s.set_uniform_color(s.c1o_loc, c1);
		s.set_uniform_color(s.c2i_loc, c2);
		s.set_uniform_color(s.c2o_loc, c2);
		s.set_uniform_color(s.c3i_loc, c3);
		s.set_uniform_color(s.c3o_loc, c3);
		s.set_uniform_float(s.rad_loc, draw_radius);
		s.set_uniform_float(s.off_loc, ((is_dynamic ? 0.0 : 100.0*pos.x) + 0.001*tfticks)); // used as a hash if static
		s.set_uniform_vector3d(s.vd_loc, (get_camera_pos() - pos).get_norm()); // local object space
		fgPushMatrix();
		translate_to(pos);
		draw_quads();
		fgPopMatrix();

		if (use_scale > 0.5) {
			shader_t cs;
			cs.set_vert_shader("no_lighting_tex_coord");
			cs.set_frag_shader("alpha_mask");
			cs.begin_shader();
			cs.add_uniform_float("min_alpha", (1.5 - use_scale));
			cs.add_uniform_int("tex0", 0);
			cs.set_cur_color(colorRGBA(0.6, 1.0, 1.0, 0.25));
			select_texture(NOISE_GEN_TEX);
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK);
			draw_subdiv_sphere(pos, draw_radius, N_SPHERE_DIV, 8, 1); // tscale=8
			glDisable(GL_CULL_FACE);
			s.enable();
		}
		// Note: tid may not be allocated if teleporter is visible in the view from the destination point
		if (tid && do_portal_draw()) { // transparent, and not too far away
			assert(tid); // must have been setup in create_portal_texture()
			bind_2d_texture(tid);
			quad_batch_draw qbd;
			qbd.add_billboard(pos, camera_pdu.pos, plus_z, WHITE, radius, radius, tex_range_t(1.0, 0.0, 0.0, 1.0)); // swap X tex coords (something backwards?)
			shader_t ts;
			ts.set_vert_shader("no_lighting_tex_coord");
			ts.set_frag_shader("circular_portal");
			ts.begin_shader();
			ts.add_uniform_float("min_alpha", 0.0);
			ts.add_uniform_int("tex0", 0);
			qbd.draw();
			ts.end_shader(); // not actually needed, but added for clarity
			s.enable(); // back to the main shader
		}
	}
	if (0 && camera_pdu.sphere_visible_test(dest, 0.25*radius)) { // draw dest (debugging)
		draw_single_colored_sphere(dest, 0.25*radius, N_SPHERE_DIV, BLUE);
		s.enable(); // back to the main shader
	}
}


void teleporter::write_to_cobj_file(std::ostream &out) const {
	out << "teleporter " << pos.raw_str() << " " << dest.raw_str() << " " << radius << " " << is_portal << " " << is_indoors << endl;
}


// jump pads

bool maybe_use_jump_pad(point &opos, vector3d &velocity, float oradius, int player_id) {

	for (auto i = jump_pads.begin(); i != jump_pads.end(); ++i) {
		if (i->maybe_jump(opos, velocity, oradius, player_id)) return 1;
	}
	return 0;
}

bool jump_pad::maybe_jump(point &opos, vector3d &obj_velocity, float oradius, int player_id) {

	if (!dist_less_than(pos, (opos - vector3d(0.0, 0.0, oradius)), radius+oradius)) return 0; // not close enough

	if (player_id == NO_SOURCE) { // object jump
		obj_velocity += velocity;
	}
	else { // player jump
		assert(player_id >= CAMERA_ID && player_id < num_smileys);
		player_state &ss(sstates[player_id]);
		if (ss.jump_time > 0) return 0; // already jumping, not touching the ground/jump pad
		ss.jump_time = 0.1*TICKS_PER_SECOND*velocity.z; // Note: only z velocity is used, so player can only jump straight up
	}
	gen_sound(SOUND_BOING, pos, 1.0, 1.0);
	last_used_tfticks = tfticks;
	return 1;
}


void draw_jump_pads() {

	if (jump_pads.empty()) return;
	shader_t s;
	s.begin_simple_textured_shader(); // no lighting
	for (auto i = jump_pads.begin(); i != jump_pads.end(); ++i) {i->draw(s);}
	s.end_shader();
}

void jump_pad::draw(shader_t &s) const {

	if (!camera_pdu.sphere_visible_test(pos, radius)) return;
	int ndiv(min(N_CYL_SIDES, max(4, int(200.0*get_zoom_scale()*sqrt(radius/distance_to_camera(pos))))));
	float const up_time(0.02), down_time(0.4);
	float const use_time(float(tfticks - last_used_tfticks)/TICKS_PER_SECOND);
	float z_pos(0.0);
	if (use_time < up_time) {z_pos = use_time/up_time;} // moving up (firing)
	else if (use_time < (up_time + down_time)) {z_pos = 1.0 - (use_time - up_time)/down_time;}
	float const len((0.2 + 0.8*z_pos)*radius), radius2(0.6*radius), thickness(0.2*radius);
	vector3d const dir(velocity.get_norm());
	s.set_cur_color(WHITE);
	select_texture(HAZARD_TEX);
	draw_fast_cylinder((pos + len*dir), (pos + (len + thickness)*dir), radius, radius, ndiv, 1, 1, 0, nullptr, 0.15);
	s.set_cur_color(GRAY);
	select_no_texture();
	draw_fast_cylinder(pos, (pos + len*dir), radius2, radius2, ndiv, 0, 0);
}

