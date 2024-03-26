// 3D World - Game Mode Weapon Drawing Code
// by Frank Gennari
// 6/24/06

#include "gameplay.h"
#include "physics_objects.h"
#include "shaders.h"
#include "draw_utils.h"


int player_dodgeball_id(-1);
vector<int> weap_cobjs;
set<int> scheduled_weapons;

extern bool keep_beams, have_indir_smoke_tex, begin_motion, enable_translocator, can_do_building_action;
extern int game_mode, window_width, window_height, frame_counter, camera_coll_id, display_mode, spectate;
extern int num_smileys, left_handed, iticks, camera_view, UNLIMITED_WEAPONS, animate2, last_inventory_frame;
extern float fticks, NEAR_CLIP, FAR_CLIP;
extern double tfticks;
extern vector3d pre_ref_cview_dir;
extern point pre_ref_camera_pos;
extern obj_type object_types[];
extern obj_group obj_groups[];
extern vector<spark_t> sparks;
extern vector<beam3d> beams;
extern vector<teleporter> teleporters[3]; // static, dynamic, in-hand
extern int coll_id[];
extern blood_spot blood_spots[];
extern player_state *sstates;


void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate);
void draw_plasmaball(point const &pos0, int shooter, shader_t &shader);
void draw_translocator(point const &pos, float radius, int ndiv, int source, shader_t &shader);


void beam3d::draw(line_tquad_draw_t &drawer) const {

	if (shooter == CAMERA_ID)  return; // camera (skip for now)
	if (intensity < TOLERANCE) return; // error?
	if (!camera_pdu.line_visible_test(pts[0], pts[1])) return;
	float const mag(sqrt(intensity));
	colorRGBA c(color);
	c.alpha = CLIP_TO_01(color.alpha*mag);
	drawer.add_line_as_tris(pts[0], pts[1], 0.01*mag, 0.01*mag, c, (distant ? colorRGBA(c, 0.0) : c));
}


void draw_beams(bool clear_at_end) {

	if (beams.empty()) return;
	//timer_t timer("Draw Beams");
	static line_tquad_draw_t drawer;
	glDepthMask(GL_FALSE);
	set_additive_blend_mode();
	for (auto i = beams.begin(); i != beams.end(); ++i) {i->draw(drawer);}
	if (clear_at_end && !keep_beams) {beams.clear();}
	drawer.draw_and_clear(0.0); // noise_scale=0.0
	set_std_blend_mode();
	glDepthMask(GL_TRUE);
}


void show_blood_on_camera() {

	enable_blend();
	shader_t s;
	s.begin_simple_textured_shader(0.5);
	select_texture(BLUR_TEX);
	quad_batch_draw qbd;
	
	for (unsigned i = 0; i < NUM_BS; ++i) {
		blood_spot &bs(blood_spots[i]);
		if (bs.time <= 0) continue;
		bs.time  -= max(1, iticks); // shrink blood spot
		bs.size  *= pow(0.99f, fticks);
		bs.pos.y -= 0.00001*sqrt(bs.size)*fticks;

		if (bs.size > 0.1) { // draw it
			float const size(0.00006*bs.size);
			qbd.add_quad_dirs(bs.pos, vector3d(size, 0.0, 0.0), vector3d(0.0, size, 0.0), colorRGBA(0.7, 0.0, 0.0));
		}
		else {
			bs.size = 0.0;
			bs.time = 0;
		}
	}
	qbd.draw();
	disable_blend();
	s.end_shader();
}


point get_final_pos(point const &pos, vector3d const &dir, float radius, float scale, float &rxy, vector3d &v_trans) {

	rxy = sqrt(dir.x*dir.x + dir.y*dir.y);

	if (rxy < TOLERANCE) { // careful of divide by zero
		v_trans = vector3d(0.0, 0.0, 2*radius);
		return (pos + v_trans);
	}
	v_trans = vector3d(-radius*dir.x/rxy, -radius*dir.y/rxy, 2*radius);
	vector3d vt(v_trans*scale);
	rotate_vector3d(vector3d(-dir.y, dir.x, 0.0), safe_acosf(-dir.z), vt);
	return (pos + vt);
}


float get_bbbat_angle(float fire_val) {
	return (-30.0 + ((fire_val > 0.95) ? 1120.0*(fire_val - 0.95) : 70.0*(0.95 - fire_val)));
}


void add_weapon_cobj(point const &pos, vector3d const &dir, float cradius, float dpos, float fire_val, int wid, int wmode) {

	if (wid == W_UNARMED) return;
	if (dir.mag() < TOLERANCE) {cout << TXT(dir.str()) << TXT(pos.str()) << TXT(wid) << TXT(wmode) << endl;}
	assert(dir != zero_vector);
	bool const DRAW_WEAP_COBJ(0); // for debugging
	int const surfs((wid == W_BLADE || wid == W_M16 || wid == W_SHOTGUN || wid == W_LASER) ? 1 : 0); // no cylinder ends
	cobj_params cp(0.8, BLACK, DRAW_WEAP_COBJ, 1, NULL, 0, -1, 1.0, surfs, 0.0, 0.0, 1); // special mode - shadow but no coll
	float rxy, radius;
	vector3d v_trans;
	point const pos0(get_final_pos(pos, -dir, cradius, 1.0, rxy, v_trans));
	bool const has_xy_comp(dir.x != 0.0 || dir.y != 0.0);

	// Note: a few weapon's cobj(s) depend on fire_val, but we could add more
	switch (wid) {
	case W_UNARMED:
		break;

	case W_BALL:
	case W_SBALL:
	case W_LANDMINE:
	case W_GRENADE:
	case W_CGRENADE:
	case W_STAR5:
	case W_SAWBLADE:
	case W_XLOCATOR: // close enough
		{
			int const oid(weapons[wid].obj_id);
			radius = 0.4*object_types[oid].radius;
			if (wid == W_CGRENADE || (wid == W_GRENADE && (wmode & 1))) radius *= 1.2;
			weap_cobjs.push_back(add_coll_sphere(pos0, radius, cp));
		}
		break;

	case W_BLADE: {
			radius = 0.032;
			float const dist(max(dpos-radius, radius)); // close enough
			assert(dist > TOLERANCE);
			point const pos1(pos0 - dir*cradius + point(0.0, 0.0, 0.01)), pos2(pos1 + dir*dist);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.0025, 0.0025, cp));
			vector3d vrot(dir);
			if (has_xy_comp) rotate_vector3d(cross_product(dir, plus_z), PI_TWO, vrot);
			if (fire_val > 0.0) rotate_vector3d(dir, -540.0*fire_val/TO_DEG, vrot);
			assert(vrot.mag() > TOLERANCE);
			cp.surfs = 0;
			weap_cobjs.push_back(add_coll_cylinder(pos2, (pos2 + vrot*(0.01*radius)), radius, radius, cp));
		}
		break;

	case W_ROCKET:
	case W_SEEK_D:
	case W_RAPTOR:
		{
			radius = 0.95*object_types[weapons[wid].obj_id].radius;
			float const rscale((wid == W_SEEK_D) ? 4.8 : ((wid == W_RAPTOR) ? 8.8 : 5.8)), cylin_radius(((wid == W_RAPTOR) ? 0.9 : 0.8)*radius);
			point const pos1(pos0 + dir*radius), pos2(pos0 - dir*(rscale*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, cylin_radius, cylin_radius, cp));
		}
		break;

	case W_PLASMA: { // Note: plasma ball itself is not a collision object
			point const pos1(pos0 - dir*0.15), pos2(pos0 + dir*0.03);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.01,  0.01, cp));
			weap_cobjs.push_back(add_coll_cylinder(pos2, (pos1 + dir*0.25), 0.005, 0.0,  cp));
		}
		break;

	case W_M16: {
			point pos1(pos + vector3d(pos0, pos)*0.6), pos2(pos1);
			pos1 -= dir*0.072;

			if ((wmode&1) == 0) { // normal
				radius = 0.0025;
				pos2  += dir*0.076;
				weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 2.8*radius, 1.5*radius, cp));
			}
			else { // shrapnel chaingun
				radius = 0.0048;
				pos2  += dir*0.04;
				weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 2.0*radius, 2.0*radius, cp));
			}
		}
		break;

	case W_SHOTGUN: {
			radius = 0.0042;
			point const pos1(pos + vector3d(pos0, pos)*0.6 - dir*0.072), pos2(pos1 + dir*0.121);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 2.0*radius, 2.0*radius, cp));
		}
		break;

	case W_BBBAT: { // this was a real bitch
			radius = 0.004;
			vector3d dir2(dir);
			if (has_xy_comp) rotate_vector3d(vector3d(-dir.y, dir.x, 0.0), PI/4.0, dir2);
			point const pos1(pos + dir*(SQRT2*cradius) - dir2*cradius);
			vector3d vr((left_handed ? 0.5 : -0.5), 0.5, 0.0);
			rotate_vector3d(vector3d(0.0, 0.0, -1.0), atan2(dir.y, dir.x), vr);
			rotate_norm_vector3d_into_plus_z(dir2, vr, -1.0); // inverse rotate
			rotate_vector3d(vr, get_bbbat_angle(fire_val)/TO_DEG, dir2);
			assert(dir2.mag() > TOLERANCE);
			point const pos2(pos1 + dir2*(16.0*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, radius, 1.2*radius, cp));
			weap_cobjs.push_back(add_coll_sphere(pos2, 1.2*radius, cp));
		}
		break;

	case W_LASER: {
			point const pos1(pos0 - dir*0.08), pos2(pos0 + dir*0.08);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.006, 0.0015, cp));
		}
		break;

	case W_GASSER: {
			radius = 0.14*weapons[W_GASSER].blast_radius;
			point const pos1(pos0 + dir*radius), pos2(pos0 - dir*(16.0*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, radius, radius, cp));
		}
		break;

	default:
		assert(0);
	}
}


void clear_weap_cobjs() { // when should this be called?

	for (unsigned i = 0; i < weap_cobjs.size(); ++i) {
		remove_coll_object(weap_cobjs[i]);
	}
	weap_cobjs.clear();
}


point get_sstate_draw_pos(int shooter) {

	point pos(get_sstate_pos(shooter));
	
	if (shooter != CAMERA_ID && !(obj_groups[coll_id[SMILEY]].get_obj(shooter).flags & CAMERA_VIEW)) {
		pos.z += 0.8*object_types[SMILEY].radius;
	}
	return pos;
}


void player_state::update_weapon_cobjs(int i) {

	point const pos(get_sstate_draw_pos(i));
	dpos = 0.0;
		
	if (weapon == W_BLADE) {
		cb_pos = pos;
		do_cblade_damage_and_update_pos(cb_pos, i);
	}
	if (world_mode != WMODE_GROUND) return;
	if (powerup == PU_INVISIBILITY) return; // smiley's/player's shadow will still be there
	float const fire_val(float(fire_frame)/float(max(1U, weapons[weapon].fire_delay)));
	add_weapon_cobj(pos, get_sstate_dir(i), object_types[SMILEY].radius, dpos, fire_val, weapon, wmode);
	add_weapon_lights(i);
}


void update_weapon_cobjs() { // and update cblade and lighting

	clear_weap_cobjs();
	if (!game_mode) return;
	assert(sstates != NULL);
	bool const smileys_enabled(begin_motion && obj_groups[coll_id[SMILEY]].enabled);

	for (int i = CAMERA_ID; i < num_smileys; ++i) { // if invisible, don't draw the weapon
		if ((i != CAMERA_ID && !smileys_enabled) || (i == CAMERA_ID && camera_view)) continue;
		sstates[i].update_weapon_cobjs(i);
	}
}


int select_dodgeball_texture(int shooter) {

	if (UNLIMITED_WEAPONS && game_mode == GAME_MODE_DODGEBALL && !obj_groups[coll_id[BALL]].reorderable) { // can change when other players throw a ball
		if (player_dodgeball_id < 0) {player_dodgeball_id = (obj_groups[coll_id[BALL]].choose_object(1) % NUM_DB_TIDS);} // peek=1
		return dodgeball_tids[player_dodgeball_id]; // choose once and don't change - may throw a ball of a different color
	}
	bool const default_tex(sstates == NULL || shooter == NO_SOURCE || sstates[shooter].balls.empty());
	assert(shooter == NO_SOURCE || (shooter >= CAMERA_ID && shooter < num_smileys));
	return dodgeball_tids[default_tex ? 0 : (sstates[shooter].balls.back() % NUM_DB_TIDS)];
}


void rotate_into_camera_dir(point const &pos, vector3d const &dir) {
	rotate_by_vector(-dir); // undo rotation
	rotate_towards_camera(pos);
}

void draw_chaingun_section(float tx, float ty, float tz, float length, float radius, int ndiv) {
	draw_cylinder_at(point(tx, ty, tz), length, radius, radius, 2*ndiv);
	draw_circle_normal(0.0, radius, ndiv, 1, point(tx, ty, tz+0.75*length));
}
void draw_raptor_section(float tx, float ty, float tz, float length, float radius, int ndiv) {
	draw_cylinder_at(point(tx, ty, tz), length, radius, radius, 2*ndiv);
	draw_subdiv_sphere(point(tx, ty, tz       ), radius, 2*ndiv, 0, 1);
	draw_subdiv_sphere(point(tx, ty, tz+length), radius, 2*ndiv, 0, 1);
}


// Note: shader is passed all the way into here just to set the min_alpha for cblade
void draw_weapon(point const &pos, vector3d dir, float cradius, int cid, int wid, int wmode, int fframe, int p_loaded,
				 int ammo, int rot_counter, unsigned delay, int shooter, int sb_tex, float alpha, float dpos,
				 float fire_val, float scale, bool draw_pass, shader_t &shader, bool fixed_lod=0)
{
	bool const just_fired(abs(fframe - (int)delay) <= 1), is_camera(shooter == CAMERA_ID);
	
	if (draw_pass == 1) {
		if (wid == W_PLASMA && p_loaded) {}
		else if ((wid == W_M16 || wid == W_SHOTGUN) && just_fired) {}
		else return; // nothing to draw
	}
	dir.negate(); // used to be backwards
	float radius, rxy;
	//glDisable(GL_DEPTH_TEST);
	enable_blend(); // not always necessary
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(scale);
	rotate_by_vector(dir);
	assert(wid < NUM_WEAPONS);
	vector3d v_trans;
	point const pos0(get_final_pos(pos, dir, cradius, scale, rxy, v_trans));
	float const tx(-cradius*dir.x/rxy), ty(-cradius*dir.y/rxy);
	select_texture(WHITE_TEX);
	
	if (draw_pass == 0) { // draw solid objects
		int ndiv(fixed_lod ? 24 : int(get_zoom_scale()*300.0f*cradius/(distance_to_camera(pos) + SMALL_NUMBER)));
		ndiv = min(N_SPHERE_DIV, max(3, ndiv));
		unsigned const oid(weapons[wid].obj_id);
		//int const cobj(weapons[wid].need_weapon ? cid : -1);
		float rot_angle;
		bool do_texture(0);

		if (oid != UNDEF) {
			assert(oid < NUM_TOT_OBJS);
			do_texture = (object_types[oid].tid >= 0);
		}
		switch (wid) {
		case W_UNARMED:
			break;

		case W_BALL:
			if (game_mode == GAME_MODE_FPS && (wmode & 1)) { // teleporter gun is secondary fire of dodgeball
				teleporter tp;
				tp.pos    = pos0;
				tp.radius = 0.3*object_types[oid].radius;
				teleporters[2].push_back(tp); // will be drawn later
				teleporters[2].back().setup();
				break;
			}
		case W_SBALL:
		case W_LANDMINE:
		case W_GRENADE: // draw_grenade()?
		case W_CGRENADE:
			radius = 0.4*object_types[oid].radius;
			if (wid == W_CGRENADE || (wid == W_GRENADE && (wmode & 1))) {set_gold_material(shader, alpha); radius *= 1.2;} // cgrenade
			else if (wid == W_GRENADE ) {set_copper_material(shader, alpha);} // grenade
			else {shader.set_cur_color(colorRGBA(object_types[oid].color, alpha)); shader.set_specular(0.8, 40.0);}
			translate_to(v_trans);
			
			if (do_texture) {
				select_texture((wid == W_BALL) ? select_dodgeball_texture(shooter) : object_types[oid].tid);
				rotate_to_dir(dir, 90.0, 1.0); // cancel out texture rotation with camera
				fgRotate(((wid == W_BALL) ? 135.0 : 45.0), 1.0, 0.0, 0.0); // rotate the texture to face the player
			}
			if (wid == W_BALL) {draw_cube_mapped_sphere(all_zeros, radius, ndiv, do_texture);} else {draw_sphere_vbo(all_zeros, radius, ndiv, do_texture);}
			shader.clear_specular();
			break;

		case W_XLOCATOR:
			radius = 0.4*object_types[oid].radius;
			translate_to(v_trans);
			rotate_to_dir(dir, 90.0, 1.0); // cancel out rotation with camera
			fgRotate(90.0, 1.0, 0.0, 0.0); // make it vertical

			if (sstates != nullptr && (shooter == NO_SOURCE || sstates[shooter].p_ammo[W_XLOCATOR] > 0)) { // have the translocator
				draw_translocator(all_zeros, radius, N_SPHERE_DIV, shooter, shader);
			}
			else { // don't have the translocator
				shader.set_cur_color(WHITE);
				select_texture(HAZARD_TEX);
				draw_cube(all_zeros, 1.6*radius, 1.0*radius, 2.4*radius, 1, 0.5); // hazard textured remote control box
				select_texture(WHITE_TEX);
				shader.set_cur_color(colorRGBA(1.0, 0.15, 0.0, 1.0)); // reddish orange
				draw_sphere_vbo(vector3d(0.0, -0.4*radius, 0.5*radius), 0.48*radius, 32, 0); // fire button
				shader.set_cur_color(RED);
				shader.set_color_e(RED);
				draw_sphere_vbo(vector3d(0.5*radius, -0.5*radius, 1.05*radius), 0.1*radius, 24, 0); // power light
				shader.clear_color_e();
				set_silver_material(shader, alpha);
				fgTranslate(0.0, 0.0, 1.4*radius);
				draw_cylinder(2.0*radius, 0.12*radius, 0.12*radius, 24, 0); // antenna
				draw_sphere_vbo(vector3d(0.0, 0.0, 2.0*radius), 0.2*radius, 24, 0);
				shader.clear_specular();
			}
			break;

		case W_STAR5:
			radius = object_types[oid].radius;
			translate_to(v_trans);
			rotate_to_dir(dir, 90.0, 1.0);  // cancel out texture rotation with camera
			fgRotate(45.0, 1.0, 0.0, 0.0); // rotate the texture to face the player
			shader.set_cur_color(colorRGBA(object_types[oid].color, alpha));
			shader.set_specular(0.8, 40.0);
			draw_star(all_zeros, plus_z, zero_vector, radius, 0.0, 0); // Note: +z may not be the correct normal?
			shader.clear_specular();
			break;

		case W_BLADE:
		case W_SAWBLADE:
			{
				static int lfc(0);
				static float angle(0.0);
				
				if (frame_counter != lfc && animate2) {
					lfc    = frame_counter;
					angle += 6.0*fticks;
				}
				radius = 0.032;
				translate_to(v_trans);
				rotate_to_dir(dir, 90.0, 1.0); // cancel out texture rotation with camera
				fgPushMatrix();
				fgTranslate(0.0, -0.009, -0.025);
				fgScale(1.0, 1.0, -1.0);
				int const ndiv2(max(4, 2*(is_camera ? N_SPHERE_DIV : ndiv)/3));
				float const len(dpos + radius), dv(4.0*radius/ndiv2);
				assert(dv > 0.0);
				shader.set_cur_color(colorRGBA(0.49, 0.51, 0.53, alpha));
				draw_fast_cylinder(point(0.0, 0.0, len), all_zeros, 0.0025, 0.0025, ndiv2, 0, 0);
				fgPopMatrix();
				fgRotate(90.0, 1.0, 0.0, 0.0); // into horizontal plane
				fgTranslate(0.0, -0.028, 0.01);
				if (fframe > 0 && !(wmode&1)) fgRotate((540.0*fframe)/delay, 0.0, 1.0, 0.0);
				fgRotate(angle, 0.0, 0.0, 1.0);
				shader.set_cur_color(colorRGBA(WHITE, alpha));
				shader.add_uniform_float("min_alpha", 0.95*alpha);
				shader.set_specular(0.9, 90.0);
				float dz((ammo > 1) ? -0.025*radius*ammo : 0.0);
				select_texture(sb_tex ? (int)SAW_B_TEX : (int)SAW_TEX);

				for (int w = 0; w < max(1, ammo); ++w) { // draw a blade for each ammo
					draw_tquad(radius, radius, dz);
					dz += 0.05*radius;
				}
				shader.clear_specular();
				shader.add_uniform_float("min_alpha", 0.01);
			}
			break;

		case W_ROCKET:
			radius = 0.95*object_types[ROCKET].radius;
			select_texture(get_texture_by_name("metal_plate.jpg"));
			set_copper_material(shader, alpha, 2.0); // bright copper
			rot_angle = max(0.0f, 10.0f*(fire_val - 0.7f)); // recoil
			fgRotate(rot_angle, -dir.y, dir.x, 0.0); // could probably use rotate_into_plus_z()
			fgTranslate(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0);
			draw_cylinder(6.8*radius, 0.8*radius, 0.8*radius, 2*ndiv);
			draw_circle_normal(0.0, 0.8*radius, ndiv, 1, 5.0*radius);
			select_texture(WHITE_TEX);
			// draw the sight
			fgTranslate(0.8*radius, 0.0, 6.5*radius);
			fgRotate(90.0, 0.0, 1.0, 0.0);
			if (wmode&1) {shader.set_cur_color(colorRGBA(BLACK, alpha));} else {set_gold_material(shader, alpha);} // black/gold
			draw_cylinder(0.4*radius, 0.15*radius, 0.0, ndiv);
			shader.clear_specular();
			break;

		case W_SEEK_D: // similar to rocket
			radius = 0.95*object_types[SEEK_D].radius;
			select_texture(SHIP_HULL_TEX);
			set_brass_material(shader, alpha, 2.0); // bright brass
			rot_angle = max(0.0f, 15.0f*(fire_val - 0.8f)); // recoil
			fgRotate(rot_angle, -dir.y, dir.x, 0.0);
			fgTranslate(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0);
			draw_cylinder(5.8*radius, 0.8*radius, 0.8*radius, 2*ndiv);
			draw_circle_normal(0.0, 0.8*radius, ndiv, 1, 4.0*radius);
			shader.clear_specular();
			select_texture(WHITE_TEX);
			break;

		case W_RAPTOR: { // similar to rocket
			radius = 0.95*object_types[(wmode&1) ? (unsigned)FREEZE_BOMB : (unsigned)RAPT_PROJ].radius;
			if (wmode&1) {shader.set_cur_color(colorRGBA(FREEZE_COLOR, alpha)); shader.set_specular(0.8, 80.0);} // freeze mode
			else {set_gold_material(shader, alpha);}
			rot_angle = max(0.0f, 8.0f*(fire_val - 0.6f)); // recoil
			fgRotate(rot_angle, -dir.y, dir.x, 0.0);
			fgTranslate(1.0*tx, 1.0*ty, 0.0);
			draw_cylinder_at(all_zeros, 8.8*radius, 0.9*radius, 0.9*radius, 2*ndiv);
			//draw_circle_normal(0.0, 0.8*radius, ndiv, 1, point(0.0, 0.0, 7.0*radius));
			draw_cylinder_at(point(0.0, 0.0, 8.8*radius), -1.8*radius, 0.9*radius, 0.0, 2*ndiv); // cone
			// draw shell chambers
			radius *= 0.5;
			float const rdx(1.5*radius*dir.x/rxy), rdy(1.5*radius*dir.y/rxy), length(5.0*radius), tz(9.0*radius);
			if (wmode&1) {set_silver_material(shader, alpha);} else {set_copper_material(shader, alpha);}
			fgPushMatrix();
			fgRotate(15.0*rot_counter, 0.0, 0.0, 1.0);
			draw_raptor_section( rdx,  rdy, tz, length, radius, ndiv);
			draw_raptor_section(-rdx, -rdy, tz, length, radius, ndiv);
			fgRotate(90.0, 0.0, 0.0, 1.0);
			draw_raptor_section( rdx,  rdy, tz, length, radius, ndiv);
			draw_raptor_section(-rdx, -rdy, tz, length, radius, ndiv);
			fgPopMatrix();
			shader.clear_specular();
			break;
		}
		case W_PLASMA:
			shader.set_cur_color(colorRGBA(BLACK, alpha));
			shader.set_specular(0.8, 10.0);
			rot_angle = max(0.0f, 2.0f*(fire_val - 0.7f));
			fgRotate(rot_angle, -dir.y, dir.x, 0.0);
			fgTranslate(0.0, 0.0, 0.15);

			if (shooter != NO_SOURCE) {
				draw_circle_normal(0.017, 0.0175, 2*ndiv, 1);
				rotate_to_dir(dir, 0.0, 1.0);
				draw_fast_cylinder(point(-0.014, 0.01, 0.0), point(-0.01, -0.001, -0.15), 0.0001, 0.0001, ndiv/2, 0);
				rotate_to_dir(dir, 0.0, -1.0);
			}
			shader.set_cur_color(colorRGBA(RED, alpha));
			draw_circle_normal(0.004, 0.0043, ndiv, 1);
			fgTranslate(tx, ty, -0.15);
			set_gold_material(shader, alpha);
			draw_circle_normal(0.0075, 0.009, ndiv, 1, 0.15);
			shader.set_cur_color(colorRGBA(BLACK, alpha));
			shader.set_specular(0.8, 10.0);
			draw_cylinder(0.15, 0.005, 0.005, 2*ndiv);
			fgPushMatrix();
			fgTranslate(0.0, 0.0, 0.15);
			draw_cylinder(0.07, 0.005, 0.0, 2*ndiv);
			//draw_cylinder(0.18, radius, radius, 2*ndiv);
			{
				colorRGBA color(0.8, 0.6, 1.0, 0.5*alpha);
				if (p_loaded) {shader.set_black_diffuse_emissive_color(color);} else {shader.set_cur_color(color);} // emissive when loaded
				fgScale(1.0, 1.0, 0.2);
				for (unsigned i = 0; i < 3; ++i) {draw_sphere_vbo_back_to_front(point(0.0, 0.0, -0.25+0.08*i), 0.01, ndiv, 0);}
				if (p_loaded) {shader.clear_color_e();}
			}
			fgPopMatrix();
			shader.clear_specular();
			break;

		case W_M16:
			if ((wmode&1) == 0) { // normal
				radius = 0.0025;
				shader.set_cur_color(colorRGBA(0.04, 0.04, 0.04, alpha));
				shader.set_specular(0.8, 50.0);
				rot_angle = max(0.0, 1.0*fire_val);
				fgRotate(rot_angle, -dir.y, dir.x, 0.0);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.076), 0.064,     radius,     radius, 2*ndiv, 1);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.000), 0.076, 2.8*radius, 2.0*radius, 2*ndiv, 1);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.136), 0.012, 1.5*radius, 1.5*radius, 2*ndiv, 1);
			}
			else { // shrapnel chaingun
				radius = 0.004;
				float const rdx(1.4*radius*dir.x/rxy), rdy(1.4*radius*dir.y/rxy), length(0.11);
				set_silver_material(shader, alpha);
				fgTranslate(0.6*tx, 0.6*ty, 0.0);
				fgPushMatrix();
				fgRotate(15.0*rot_counter, 0.0, 0.0, 1.0);
				draw_chaingun_section( rdx,  rdy, 0.0, length, radius, ndiv);
				draw_chaingun_section(-rdx, -rdy, 0.0, length, radius, ndiv);
				fgRotate(90.0, 0.0, 0.0, 1.0);
				draw_chaingun_section( rdx,  rdy, 0.0, length, radius, ndiv);
				draw_chaingun_section(-rdx, -rdy, 0.0, length, radius, ndiv);
				shader.set_cur_color(colorRGBA(0.3, 0.3, 0.3, alpha));
				draw_cylinder_at(point(0.0, 0.0, 0.08), 0.004, 2.45*radius, 2.45*radius, 2*ndiv, 1); // outer band
				fgPopMatrix();
			}
			shader.clear_specular();
			break;

		case W_SHOTGUN:
			{
				radius = 0.0045;
				float const rdx(radius*dir.x/rxy), rdy(radius*dir.y/rxy);
				set_gold_material(shader, alpha);
				rot_angle = max(0.0, 8.0*fire_val);
				fgRotate(rot_angle, -dir.y, dir.x, 0.0);
				fgTranslate(0.6*tx, 0.6*ty, 0.0);
				fgRotate(90.0, 0.0, 0.0, 1.0);
				point const translates[2] = {point(rdx, rdy, 0.0), point(-0.9*rdx, -0.9*rdy, 0.0)};
				
				for (unsigned i = 0; i < 2; ++i) {
					draw_cylinder_at(translates[i], 0.12, radius, radius, 2*ndiv);
					draw_circle_normal(0.0, radius, ndiv, 1, translates[i]+vector3d(0.0, 0.0, 0.09));
				}
				shader.clear_specular();
			}
			break;
			
		case W_BBBAT:
			radius = 0.004;
			shader.set_cur_color(colorRGBA(LT_BROWN, alpha));
			select_texture(is_camera ? (int)PLAYER_BBB_TEX : (int)WOOD_TEX); // customize the player's baseball bat
			fgRotate(45.0, -dir.y, dir.x, 0.0);
			fgTranslate(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0); // cancel out texture rotation with camera
			fgRotate(get_bbbat_angle(fire_val), (left_handed ? 0.5 : -0.5), 0.5, 0.0);
			draw_cylinder(16.0*radius, radius, 1.2*radius, ndiv);
			draw_sphere_vbo(point(0.0, 0.0, 16.0*radius), 1.2*radius, ndiv, 1);
			break;

		case W_LASER:
			if (shooter == CAMERA_ID && fire_val > 0.0) {
				colorRGBA const laser_color(get_laser_beam_color(shooter));
				//beams.push_back(beam3d(0, shooter, (pos0 + dir*0.08), (pos0 + dir*1.08), laser_color)); // should probably use this instead
				fgPushMatrix();
				shader.set_black_diffuse_emissive_color(laser_color);
				fgTranslate(0.0, 0.0, 0.148);
				fgRotate(30.0, -dir.y, dir.x, 0.0);
				draw_cylinder_at(point(tx, ty, 0.0), 0.103, 0.0007, 0.0007, ndiv);
				shader.clear_color_e();
				fgPopMatrix();
			}
			shader.set_cur_color(colorRGBA(0.45, 0.45, 0.45, alpha));
			shader.set_specular(0.8, 50.0);
			draw_cylinder_at(point(tx, ty, 0.04), 0.16, 0.006, 0.0015, 2*ndiv);
			draw_sphere_vbo(point(tx, ty, 0.04), 0.006, ndiv, 0);
			shader.clear_specular();
			break;

		case W_GASSER:
			radius = 0.14*weapons[W_GASSER].blast_radius;
			shader.set_cur_color(colorRGBA(OLIVE*0.7, alpha));
			shader.set_specular(0.7, 30.0);
			draw_cylinder_at(point(tx, ty, 0.0), 16.0*radius, radius, radius, 2*ndiv);
			draw_circle_normal(0.0, radius, ndiv, 1, point(tx, ty, 8.0*radius));
			shader.clear_specular();
			break;

		default:
			cout << "Error: No draw model for weapon: " << wid << endl;
			assert(0);
		}
	}
	else { // draw alpha blended objects
		switch (wid) {
		case W_PLASMA:
			if (p_loaded && shooter != NO_SOURCE) {
				fgTranslate(tx, ty, 0.0);
				rotate_to_dir(dir, 0.0, 1.0); // cancel out texture rotation with camera
				draw_plasmaball(pos, shooter, shader);
			}
			break;

		case W_M16:
			if (just_fired) {
				float const size(((wmode&1) == 0) ? 0.02 : 0.0272);
				shader.set_black_diffuse_emissive_color(ORANGE);
				fgTranslate(0.6*tx, 0.6*ty, (((wmode&1) == 0) ? 0.15 : 0.12));
				if (!is_camera) rotate_into_camera_dir(pos0, dir); // pos0 is approximate
				set_additive_blend_mode();
				select_texture(FLARE1_TEX);
				set_std_blend_mode();
				draw_tquad(size, size, 0.0);
				shader.clear_color_e();
			}
			break;

		case W_SHOTGUN:
			if (just_fired) {
				radius = 0.0042;
				float const rdx(radius*dir.x/rxy), rdy(radius*dir.y/rxy);
				shader.set_black_diffuse_emissive_color(ORANGE);
				set_additive_blend_mode();
				fgTranslate(0.6*tx, 0.6*ty, 0.0);
				fgRotate(90.0, 0.0, 0.0, 1.0);
				point const translates[2] = {point(-0.9*rdx, -0.9*rdy, 0.124), point(1.9*rdx, 1.9*rdy, -0.002)};
				select_texture(FLARE2_TEX);
				
				for (unsigned i = 0; i < 2; ++i) {
					translate_to(translates[i]);
					fgPushMatrix();
					if (!is_camera) {rotate_into_camera_dir(pos0, dir);} // pos0 is approximate
					draw_tquad(8.0*radius, 8.0*radius, 0.0); // can't rotate towards camera, already rotated
					fgPopMatrix();
				}
				set_std_blend_mode();
				shader.clear_color_e();
			}
			break;
		}
	}
	select_texture(WHITE_TEX);
	fgPopMatrix();
	disable_blend();
	//glEnable(GL_DEPTH_TEST);
}


void draw_weapon_simple(point const &pos, vector3d const &dir, float radius, int cid, int wid, float scale, shader_t &shader, int shooter, bool fixed_lod, float apha) {
	draw_weapon(pos, dir, radius, cid, wid, 0, 0, 0, 1, 0, 2, shooter, 0, apha, 0.0, 0.0, scale, 0, shader, fixed_lod);
}


bool weap_has_transparent(int shooter) {

	assert(shooter == CAMERA_ID || shooter < num_smileys);
	if (!game_mode)      return 0;
	if (sstates == NULL) return 0; // not initialized - should this be an error?
	player_state &sstate(sstates[shooter]);
	if (sstate.powerup == PU_INVISIBILITY) return 1;

	switch (sstate.weapon) {
	case W_PLASMA:
		return (1 || sstate.plasma_loaded); // charging rings are semi-transparent
	case W_M16:
	case W_SHOTGUN:
		return (abs(sstate.fire_frame - (int)weapons[sstate.weapon].fire_delay) <= 1);
	}
	return 0;
}


int get_shooter_coll_id(int shooter) {
	return ((shooter == CAMERA_ID) ? camera_coll_id : obj_groups[coll_id[SMILEY]].get_obj(shooter).coll_id);
}


void draw_weapon_in_hand_real(int shooter, bool draw_pass, shader_t &shader, int reflection_pass=0) {

	assert(shooter == CAMERA_ID || shooter < num_smileys);
	assert(sstates != NULL);
	player_state &sstate(sstates[shooter]);
	int const wid(sstate.weapon);
	if (wid == W_UNARMED) return;
	float alpha(1.0);
	float const cradius(object_types[SMILEY].radius);
	vector3d dir(get_sstate_dir(shooter));
	bool cull_face(0);

	if (shooter == CAMERA_ID) { // camera/player
		if (camera_view) return;

		if (sstate.powerup == PU_INVISIBILITY) {
			bool const flash(sstate.powerup_time < int(4*TICKS_PER_SECOND) && ((4*sstate.powerup_time/TICKS_PER_SECOND)&1));
			alpha *= (flash ? 0.8 : 0.4);

			if (wid != W_BLADE) {
				cull_face = 1;
				glEnable(GL_CULL_FACE);
				glCullFace(GL_BACK);
			}
		}
		if (reflection_pass) {dir = pre_ref_cview_dir;} // use player dir, not reflection camera dir
		else if (0) { // add weapon sway (disabled)
			static vector3d prev_dir;
			float const delta_dir(2.5*(1.0 - pow(0.7f, fticks)));
			dir = delta_dir*dir + (1.0 - delta_dir)*prev_dir;
			if (draw_pass == 1) {prev_dir = dir;} // update on second draw pass so that it's consistent across Z-prepass, transparent, and opaque passes
		}
	}
	else { // smiley - out of sync by a frame?
		assert(shooter >= 0 && shooter < num_smileys);
		if (sstate.powerup == PU_INVISIBILITY) return;
		if (sstate.weapon == W_BALL && game_mode == GAME_MODE_DODGEBALL) return; // dodgeball already drawn
		reflection_pass = 0; // irrelevant for smileys
	}
	int const cid(get_shooter_coll_id(shooter));
	unsigned const delay(max(1u, weapons[wid].fire_delay));
	float const fire_val((float)sstate.fire_frame/(float)delay);
	point const pos((draw_pass == 0 && wid == W_BLADE) ? sstate.cb_pos : (reflection_pass ? pre_ref_camera_pos : get_sstate_draw_pos(shooter)));
	float const scale(CAMERA_RADIUS/0.06); // weapons scale with camera radius (Note: inverse scale is pre-applied to cradius so that it cancels out inside draw_weapon())
	select_texture(WHITE_TEX); // always textured
	draw_weapon(pos, dir, cradius/scale, cid, wid, sstate.wmode, sstate.fire_frame, sstate.plasma_loaded, sstate.p_ammo[wid],
		sstate.rot_counter, delay, shooter, (sstate.cb_hurt > 20), alpha, sstate.dpos, fire_val, scale, draw_pass, shader);
	if (cull_face) {glDisable(GL_CULL_FACE);}
}


void draw_weapon_in_hand(int shooter, shader_t &shader, int reflection_pass) {

	if (!game_mode) return;
	draw_weapon_in_hand_real(shooter, 0, shader, reflection_pass);
	if (reflection_pass) {draw_weapon_in_hand_real(shooter, 1, shader, reflection_pass);} // draw the other pass now
	else {scheduled_weapons.insert(shooter);} // should not be duplicates, but just in case draw_scheduled_weapons() isn't called
}

void draw_camera_weapon(bool want_has_trans, int reflection_pass) {

	if (!game_mode || (weap_has_transparent(CAMERA_ID) != want_has_trans) || (game_mode == GAME_MODE_DODGEBALL && world_mode != WMODE_GROUND)) return;
	if (reflection_pass && !camera_pdu.sphere_visible_test(pre_ref_camera_pos, 2.0*CAMERA_RADIUS)) return; // player + weapon not visible in reflection
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0.0, 0.0, 0, 0, 0, 0); // no rain/snow
	draw_weapon_in_hand(CAMERA_ID, s, reflection_pass);
	s.end_shader();

	if (!want_has_trans && !reflection_pass) {
		// to make the weapon always in front, draw it again with a custom depth value just in front of the near plane;
		// doesn't work for plasma cannon due to alpha blend; alpha test is needed for cblade
		s.set_vert_shader("pos_only");
		s.set_frag_shader("write_const_depth");
		s.begin_shader();
		s.add_uniform_float("depth_value", 1.1*NEAR_CLIP);
		draw_weapon_in_hand(CAMERA_ID, s, 0);
		s.end_shader();
	}
}


void draw_scheduled_weapons(bool clear_after_draw) {
	
	if (scheduled_weapons.empty()) return;
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0.0, 0.0, 0, 0, 0, 0); // keep_alpha=1; no rain/snow enable smoke?
	for (auto i = scheduled_weapons.begin(); i != scheduled_weapons.end(); ++i) {draw_weapon_in_hand_real(*i, 1, s);}
	s.end_shader();
	if (clear_after_draw) {scheduled_weapons.clear();}
}


void draw_plasmaball(point const &pos0, int shooter, shader_t &shader) { // and shoot lightning

	if (shooter == NO_SOURCE) return;
	assert(sstates != NULL);
	int min_i(NO_SOURCE), ndiv(N_SPHERE_DIV/2);
	int const cid(coll_id[SMILEY]);
	float psize(sstates[shooter].plasma_size);
	point pos(0.0, 0.0, 0.22), spos(get_sstate_pos(shooter));
	point const player(get_camera_pos());
	enable_blend();
	glEnable(GL_CULL_FACE);
	select_texture(PLASMA_TEX);
	float radius(object_types[PLASMA].radius);
	obj_group &objg(obj_groups[cid]);
	if (shooter == CAMERA_ID || (objg.get_obj(shooter).flags & CAMERA_VIEW)) ndiv *= 3;
	draw_plasma(pos, (pos + pos0), radius, psize, ndiv, 1, 0, 0, shader);
	select_texture(WHITE_TEX);
	glDisable(GL_CULL_FACE);
	if (psize < 0.9*MAX_PLASMA_SIZE) return;

	if (is_underwater(spos, 1)) { // under water - suicide from lightning
		smiley_collision(shooter, shooter, zero_vector, pos0, PLASMA_LT_DAMAGE, PLASMA_LT_D);
	}
	// lightning eminating from plasma
	radius *= psize;
	vector3d target(zero_vector);
	float min_dist(4.0*radius + object_types[SMILEY].radius);
			
	for (unsigned i = 0; i < objg.end_id; ++i) { // test smileys
		if (objg.get_obj(i).disabled() || (int)i == shooter || same_team(i, shooter)) continue;
		float const dist(p2p_dist(spos, objg.get_obj(i).pos));

		if (dist < min_dist) { // close enough to hit a smiley?
			min_dist = dist;
			min_i    = i;
		}
	}
	if (shooter != CAMERA_ID && !same_team(CAMERA_ID, shooter)) { // test player
		if (p2p_dist(spos, player) < min_dist) {min_i = CAMERA_ID;} // close enough to hit a the player?
	}
	if (min_i >= CAMERA_ID) {
		vector3d const v(get_sstate_dir(shooter));
		point const apos(get_sstate_pos(shooter));
		target = apos - spos;
		UNROLL_3X(target[i_] -= pos[i_]*v[i_];)
		target.normalize();
	}
	set_additive_blend_mode();
	int const num_rays((rand()%6)-2);

	for (int i = 0; i < num_rays; ++i) {
		bool const hit(min_i >= CAMERA_ID && rand_float() < 0.6);
		point pos2(pos);

		if (hit) { // targeted to a smiley or the player
			UNROLL_3X(pos2[i_] += rand_uniform(0.9, 1.3)*radius*target[i_];)
			vadd_rand(pos2, 0.5*radius);
			smiley_collision(min_i, shooter, zero_vector, pos2, PLASMA_LT_DAMAGE, PLASMA_LT_D);
		}
		else {
			vadd_rand(pos2, 1.7*radius);
		}
		unsigned const npts(rand() & 15);
		line3d line;
		line.color = LITN_C;
		line.width = (hit ? 0.3 : 0.15);
		line.points.resize(npts + 2);
		line.points[0] = pos;
		line.points[1] = pos2;

		for (unsigned j = 0; j < npts; ++j) {
			UNROLL_3X(pos2[i_] *= rand_uniform(1.01, 1.2);)
			line.points[j+2] = pos2;
		}
		line.draw_lines(1, (i+1 < num_rays)); // uses a custom shader; fade_ends=1, end draw on last ray
	}
	set_std_blend_mode();
	shader.enable(); // re-enable after line drawing (which uses a different shader)
}


void add_weapon_lights(int shooter) {

	if (sstates == NULL) return;
	assert(shooter >= CAMERA_ID && shooter < num_smileys);
	player_state const &sstate(sstates[shooter]);

	switch (sstate.weapon) {
	case W_PLASMA:
		if (sstate.plasma_loaded) {
			point pos(get_sstate_draw_pos(shooter));
			pos.z += 0.22;
			add_dynamic_light(min(3.5, 45.0*sstate.plasma_size*object_types[PLASMA].radius), pos, get_plasma_color(sstate.plasma_size));
		}
		break;
	case W_GASSER:
		if (sstate.wmode & 1) {
			add_dynamic_light(0.2, (get_sstate_draw_pos(shooter) + get_sstate_dir(shooter)*(16*0.14*weapons[W_GASSER].blast_radius)), ORANGE); // add sparks?
		}
	}
}


void show_crosshair(colorRGBA const &color, int in_zoom) {

	float const scale((world_mode == WMODE_UNIVERSE) ? 0.25 : 1.0); // closer near clip for planets
	float const xy_scale(scale*((can_do_building_action && !in_zoom) ? 1.5 : 1.0));
	float const xy1(0.0006*xy_scale), xy2(0.0002*xy_scale), zval(-0.05*scale);
	float const xy[8] = {-xy1, -xy2, xy1, xy2, 0.0, 0.0, 0.0, 0.0};
	glDisable(GL_DEPTH_TEST);

	if (in_zoom) {
		shader_t s;
		s.begin_color_only_shader(color);
		fgPushMatrix();
		fgScale(2.0, 2.0, 1.0);
		vert_wrap_t verts[8];
		for (unsigned i = 0; i < 8; ++i) {verts[i] = point(xy[i], xy[(i+4)&7], zval);}
		draw_verts(verts, 8, GL_LINES);
		fgPopMatrix();
		s.end_shader();
	}
	else {
		point_sprite_drawer psd;
		for (unsigned i = 0; i < 4; ++i) {psd.add_pt(vert_color(point(xy[2*i], xy[(2*i+4)&7], zval), color));}
		psd.add_pt(vert_color(point(0.0, 0.0, zval), color));
		psd.draw(WHITE_TEX, 1.0); // draw with points of size 1 pixel, not blended
	}
	glEnable(GL_DEPTH_TEST);
}


void draw_inventory() {

	// Note: currently only weapons as there are no other collectible pickup items
	if (!game_mode || sstates == NULL) return;
	if (spectate) return; // onlys show if player is alive and playing
	float const elapsed_secs(float(frame_counter - last_inventory_frame)/TICKS_PER_SECOND);
	if (elapsed_secs > 3.0) return; // inventory only shows up for 2-3s after item pickup or inventory change
	float alpha(1.0);
	if (elapsed_secs > 2.0) {alpha = 3.0 - elapsed_secs;} // fade over the last second
	player_state const &sstate(sstates[CAMERA_ID]);
	vector<unsigned> weapons;

	for (unsigned i = 1; i < NUM_WEAPONS; ++i) { // skip unarmed
		// if player doesn't have this weapon, or it's out of ammo, don't show it
		if (game_mode == GAME_MODE_DODGEBALL && i != W_BALL) continue; // dodgeballs only
		if ((enable_translocator && i == W_XLOCATOR) || (!sstate.no_weap_id(i) && !sstate.no_ammo_id(i))) {weapons.push_back(i);}
	}
	shader_t s;
	s.begin_simple_textured_shader(); // no lighting
	float const ar(float(window_width)/float(window_height)), dx(0.006*ar), quad_sz(0.35*dx), border_sz(1.06*quad_sz), radius(0.05);
	float const x0(-0.045*ar + 0.5*(15 - int(weapons.size()))*dx); // centered on the screen
	point pos(x0, -0.045, -10.0*NEAR_CLIP);

	// draw background boxes
	quad_batch_draw qbd;
	
	for (auto w = weapons.begin(); w != weapons.end(); ++w) {
		if ((int)*w == sstate.weapon) {qbd.add_quad_dirs(pos, border_sz*plus_x, border_sz*plus_y, colorRGBA(WHITE, alpha));} // draw selection box
		qbd.add_quad_dirs(pos, quad_sz*plus_x, quad_sz*plus_y, colorRGBA(0.1, 0.1, 0.1, alpha)); // add near-black square background
		pos.x += dx;
	}
	select_texture(WHITE_TEX);
	glDisable(GL_DEPTH_TEST);
	enable_blend();
	qbd.draw();
	disable_blend();
	glEnable(GL_DEPTH_TEST);

	// draw weapons
	// UNARMED BBBAT BALL SBALL ROCKET LANDMINE SEEK_D STAR5 M16 SHOTGUN GRENADE LASER PLASMA BLADE GASSER RAPTOR XLOCATOR
	float const scale[NUM_WEAPONS] = {0.0, 2.0, 1.0, 2.0, 0.4, 2.0, 0.4, 2.0, 0.4, 0.4, 2.0, 0.4, 0.35, 1.0, 0.4, 0.4, 2.0};
	float const dist [NUM_WEAPONS] = {0.0, 1.0, 1.3, 2.6, 0.4, 2.6, 0.4, 2.6, 0.4, 0.4, 2.6, 0.6, 0.5,  1.0, 0.4, 0.4, 2.6};
	glClear(GL_DEPTH_BUFFER_BIT);
	pos.x = x0;

	for (auto w = weapons.begin(); w != weapons.end(); ++w) {
		draw_weapon_simple((pos + dist[*w]*vector3d(0.0, -0.15*radius, 0.004)), plus_y, radius, -1, *w, 0.1*scale[*w], s, CAMERA_ID, 1, alpha); // fixed_lod=1
		pos.x += dx;
	}
}

void draw_qbd_with_textured_shader(quad_batch_draw const &qbd, int tid, float min_alpha=0.0) {

	select_texture(tid);
	shader_t s;
	s.begin_simple_textured_shader(min_alpha); // no lighting
	glDisable(GL_DEPTH_TEST);
	qbd.draw();
	glEnable(GL_DEPTH_TEST);
}

void show_player_keycards() {

	if (spectate) return; // onlys show if player is alive and playing
	if (sstates == nullptr) return;
	set<unsigned> const &keycards(sstates[CAMERA_ID].keycards);
	if (keycards.empty()) return;
	float const ar(float(window_width)/float(window_height)), s(10.0*DEF_NEAR_CLIP), dx(0.06*s*ar), quad_sz_x(0.35*dx), quad_sz_y(0.7*quad_sz_x);
	point pos(0.52*s*ar, 0.52*s, -s); // top right, extending left
	quad_batch_draw qbd;

	for (auto i = keycards.begin(); i != keycards.end(); ++i) {
		qbd.add_quad_dirs(pos, quad_sz_x*plus_x, quad_sz_y*plus_y, get_keycard_color(*i));
		pos.x -= dx;
	}
	draw_qbd_with_textured_shader(qbd, KEYCARD_TEX);
}

void show_icon_image(string const &fn, float xsize, float ysize, float xpos=0.0, vector<colorRGBA> const &colors=vector<colorRGBA>()) {
	float const ar(float(window_width)/float(window_height)), s(10.0*DEF_NEAR_CLIP), quad_sz_x(0.025*s*xsize), quad_sz_y(0.025*s*ysize);
	quad_batch_draw qbd;
	point pos((0.52 - 0.05*xpos)*s*ar, 0.52*s, -s);
	qbd.add_quad_dirs(pos, quad_sz_x*plus_x, quad_sz_y*plus_y, WHITE); // top right; dx and dy are radius values
	draw_qbd_with_textured_shader(qbd, get_texture_by_name(fn, 0, 0, 0), 0.1); // wrap_mir=0 (clamp), min_alpha=0.1

	if (!colors.empty()) { // draw colors as bars in a row below the icon
		float const pitch(2.2*quad_sz_x/max(colors.size(), size_t(2))), width(0.75*pitch), hheight(0.67*quad_sz_y);
		qbd.clear();
		pos.y -= 1.25*quad_sz_y + hheight; // shift below
		pos.x -= quad_sz_x - 0.5*pitch; // shift nearly to left edge

		for (colorRGBA const &color : colors) {
			if (color.A > 0.0) {qbd.add_quad_dirs(pos, 0.5*width*plus_x, hheight*plus_y, color);} // skip transparent color
			pos.x += pitch;
		}
		shader_t s;
		s.begin_color_only_shader();
		glDisable(GL_DEPTH_TEST);
		qbd.draw();
		glEnable(GL_DEPTH_TEST);
	}
}
void show_key_icon(vector<colorRGBA> const &colors) {show_icon_image("icons/key.png",        1.0, 0.4, 0.0, colors);} // rightmost slot
void show_flashlight_icon()                         {show_icon_image("icons/flashlight.png", 1.0, 1.0, 1.0        );} // one slot  to the left
void show_pool_cue_icon  ()                         {show_icon_image("icons/pool_cue.png",   1.0, 1.0, 2.0        );} // two slots to the left

