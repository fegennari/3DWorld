// 3D World - Game Mode Weapon Drawing Code
// by Frank Gennari
// 6/24/06

#include "gameplay.h"
#include "physics_objects.h"
#include "shaders.h"
#include "draw_utils.h"


vector<int> weap_cobjs;
set<int> scheduled_weapons;

extern bool keep_beams, have_indir_smoke_tex;
extern int game_mode, window_width, window_height, frame_counter, camera_coll_id, display_mode, begin_motion;
extern int num_smileys, left_handed, iticks, camera_view, fired, UNLIMITED_WEAPONS, animate2;
extern float fticks, tfticks;
extern obj_type object_types[];
extern obj_group obj_groups[];
extern vector<spark_t> sparks;
extern vector<beam3d> beams;
extern vector<teleporter> teleporters;
extern int coll_id[];
extern blood_spot blood_spots[];
extern player_state *sstates;


void set_emissive_only(colorRGBA const &color, shader_t &shader);
void draw_star(point const &pos, vector3d const &orient, vector3d const &init_dir, float radius, float angle, int rotate);
void draw_plasmaball(point const &pos0, int shooter, shader_t &shader);


void beam3d::draw(line_tquad_draw_t &drawer) const {

	if (shooter == CAMERA_ID)  return; // camera (skip for now)
	if (intensity < TOLERANCE) return; // error?
	float const mag(sqrt(intensity));
	colorRGBA c(color);
	c.alpha = CLIP_TO_01(color.alpha*mag);
	drawer.add_line_as_tris(pts[0], pts[1], 0.01*mag, 0.01*mag, c, (distant ? ALPHA0 : c));
}


void draw_beams() {

	if (beams.empty()) return;
	line_tquad_draw_t drawer;
	for (unsigned i = 0; i < beams.size(); ++i) {beams[i].draw(drawer);}
	if (!keep_beams) {beams.clear();}
	drawer.draw();
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


float const get_bbbat_angle(float fire_val) {

	return (-30.0 + ((fire_val > 0.95) ? 1120.0*(fire_val - 0.95) : 70.0*(0.95 - fire_val)));
}


void add_weapon_cobj(point const &pos, vector3d const &dir, float cradius, float dpos, float fire_val, int wid, int wmode) {

	if (wid == W_UNARMED) return;
	assert(dir.mag() > TOLERANCE);
	bool const DRAW_WEAP_COBJ(0); // for debugging
	int const surfs((wid == W_BLADE || wid == W_M16 || wid == W_SHOTGUN || wid == W_LASER) ? 1 : 0); // no cylinder ends
	cobj_params cp(0.8, BLACK, DRAW_WEAP_COBJ, 1, NULL, 0, -1, 1.0, surfs, 0.0, 0.0, 1); // special mode - shadow but no coll
	float rxy, radius;
	vector3d v_trans;
	point const pos0(get_final_pos(pos, -dir, cradius, 1.0, rxy, v_trans));
	bool const has_xy_comp(dir.x != 0.0 || dir.y != 0.0);

	// FIXME: animate some of these based on fire_val?
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
		{
			int const oid(weapons[wid].obj_id);
			radius = 0.4*object_types[oid].radius;
			if (wid == W_CGRENADE || (wid == W_GRENADE && (wmode & 1))) radius *= 1.2;
			weap_cobjs.push_back(add_coll_sphere(pos0, radius, cp));
		}
		break;

	case W_BLADE:
		{
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
		{
			radius = 0.95*object_types[ROCKET].radius;
			point const pos1(pos0 + dir*radius), pos2(pos0 - dir*(5.8*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.8*radius, 0.8*radius, cp));
		}
		break;

	case W_SEEK_D:
		{
			radius = 0.95*object_types[SEEK_D].radius;
			point const pos1(pos0 + dir*radius), pos2(pos0 - dir*(4.8*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.8*radius, 0.8*radius, cp));
		}
		break;

	case W_PLASMA:
		{
			radius = 0.018;
			point const pos1(pos0 - dir*0.15), pos2(pos0 + dir*0.03);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.01,  0.01, cp));
			weap_cobjs.push_back(add_coll_cylinder(pos2, (pos1 + dir*0.25), 0.005, 0.0,  cp));
		}
		break;

	case W_M16:
		{
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

	case W_SHOTGUN:
		{
			radius = 0.0042;
			point const pos1(pos + vector3d(pos0, pos)*0.6 - dir*0.072), pos2(pos1 + dir*0.121);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 2.0*radius, 2.0*radius, cp));
		}
		break;

	case W_BBBAT: // this was a real bitch
		{
			radius = 0.004;
			vector3d dir2(dir);
			if (has_xy_comp) rotate_vector3d(vector3d(-dir.y, dir.x, 0.0), PI/4.0, dir2);
			point const pos1(pos + dir*(SQRT2*cradius) - dir2*cradius);
			vector3d vr((left_handed ? 0.5 : -0.5), 0.5, 0.0);
			rotate_vector3d(vector3d(0.0, 0.0, -1.0), atan2(dir.y, dir.x), vr);
			rotate_vector3d_by_vr(plus_z, dir2, vr);
			rotate_vector3d(vr, get_bbbat_angle(fire_val)/TO_DEG, dir2);
			assert(dir2.mag() > TOLERANCE);
			point const pos2(pos1 + dir2*(16.0*radius));
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, radius, 1.2*radius, cp));
			weap_cobjs.push_back(add_coll_sphere(pos2, 1.2*radius, cp));
		}
		break;

	case W_LASER:
		{
			point const pos1(pos0 - dir*0.08), pos2(pos0 + dir*0.08);
			weap_cobjs.push_back(add_coll_cylinder(pos1, pos2, 0.006, 0.0015, cp));
		}
		break;

	case W_GASSER:
		{
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


inline void rotate_to_dir(vector3d const &dir, float vadd, float vmult) {

	fgRotate(TO_DEG*vmult*atan2(dir.y, dir.x) + vadd, 0.0, 0.0, 1.0);
}


int select_dodgeball_texture(int shooter) {

	if (UNLIMITED_WEAPONS && game_mode == 2 && !obj_groups[coll_id[BALL]].reorderable) { // can change when other players throw a ball 
		return dodgeball_tids[obj_groups[coll_id[BALL]].choose_object(1) % NUM_DB_TIDS]; // peek=1
	}
	bool const default_tex(sstates == NULL || shooter == NO_SOURCE || sstates[shooter].balls.empty());
	assert(shooter >= CAMERA_ID && shooter < num_smileys);
	return dodgeball_tids[default_tex ? 0 : (sstates[shooter].balls.back() % NUM_DB_TIDS)];
}


void rotate_into_camera_dir(point const &pos, vector3d const &dir) {

	rotate_by_vector(-dir); // undo rotation
	rotate_towards_camera(pos);
}


void draw_chaingun_section(float tx, float ty, float radius, int ndiv) {

	draw_cylinder_at(point(tx, ty, 0.0), 0.11, radius, radius, 2*ndiv);
	draw_circle_normal(0.0, radius, ndiv, 1, point(tx, ty, 0.08));
}


// Note: shader is passed all the way into here just to set the min_alpha for cblade
void draw_weapon(point const &pos, vector3d dir, float cradius, int cid, int wid, int wmode, int fframe, int p_loaded,
				 int ammo, int rot_counter, unsigned delay, int shooter, int sb_tex, float alpha, float dpos,
				 float fire_val, float scale, bool draw_pass, shader_t &shader)
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
	set_fill_mode();
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
		int ndiv(int(get_zoom_scale()*280.0*cradius/(distance_to_camera(pos) + SMALL_NUMBER)));
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
		case W_SBALL:
		case W_LANDMINE:
		case W_GRENADE: // draw_grenade()?
		case W_CGRENADE:
			if (do_texture) {select_texture((wid == W_BALL) ? select_dodgeball_texture(shooter) : object_types[oid].tid);}
			radius = 0.4*object_types[oid].radius;
			if (wid == W_CGRENADE || (wid == W_GRENADE && (wmode & 1))) radius *= 1.2;
			translate_to(v_trans);
			if (do_texture) rotate_to_dir(dir, 90.0, 1.0);  // cancel out texture rotation with camera
			if (do_texture) fgRotate(45.0, 1.0, 0.0, 0.0); // rotate the texture to face the player
			shader.set_cur_color(colorRGBA(object_types[oid].color, alpha));
			shader.set_specular(0.8, 40.0);
			draw_sphere_vbo(all_zeros, radius, ndiv, do_texture);
			shader.set_specular(0.0, 0.0);
			break;

		case W_STAR5:
			radius = object_types[oid].radius;
			translate_to(v_trans);
			rotate_to_dir(dir, 90.0, 1.0);  // cancel out texture rotation with camera
			fgRotate(45.0, 1.0, 0.0, 0.0); // rotate the texture to face the player
			shader.set_cur_color(colorRGBA(object_types[oid].color, alpha));
			shader.set_specular(0.8, 40.0);
			draw_star(zero_vector, plus_z, zero_vector, radius, 0.0, 0); // Note: +z may not be the correct normal?
			shader.set_specular(0.0, 0.0);
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
				select_texture(sb_tex ? SAW_B_TEX : SAW_TEX);

				for (int w = 0; w < max(1, ammo); ++w) { // draw a blade for each ammo
					draw_tquad(radius, radius, dz);
					dz += 0.05*radius;
				}
				shader.set_specular(0.0, 0.0);
				shader.add_uniform_float("min_alpha", 0.01);
			}
			break;

		case W_ROCKET:
			radius = 0.95*object_types[ROCKET].radius;
			shader.set_cur_color(colorRGBA(0.15, 0.15, 0.15, alpha));
			shader.set_specular(0.9, 50.0);
			rot_angle = max(0.0f, 10.0f*(fire_val - 0.7f)); // recoil
			if (rot_angle != 0.0) fgRotate(rot_angle, -dir.y, dir.x, 0.0); // could probably use rotate_into_plus_z()
			fgTranslate(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0);
			draw_cylinder(6.8*radius, 0.8*radius, 0.8*radius, 2*ndiv);
			draw_circle_normal(0.0, 0.8*radius, ndiv, 1, 5.0*radius);
			// draw the sight
			fgTranslate(0.8*radius, 0.0, 6.5*radius);
			fgRotate(90.0, 0.0, 1.0, 0.0);
			shader.set_cur_color(colorRGBA((wmode&1) ? BLACK : colorRGBA(0.9, 0.65, 0.05), alpha)); // black/gold
			draw_cylinder(0.4*radius, 0.15*radius, 0.0, ndiv);
			shader.set_specular(0.0, 0.0);
			break;

		case W_SEEK_D: // similar to rocket
			radius = 0.95*object_types[SEEK_D].radius;
			shader.set_cur_color(colorRGBA(0.05, 0.05, 0.05, alpha));
			shader.set_specular(0.7, 30.0);
			rot_angle = max(0.0f, 15.0f*(fire_val - 0.8f));
			if (rot_angle != 0.0) fgRotate(rot_angle, -dir.y, dir.x, 0.0);
			draw_cylinder_at(point(tx, ty, 0.0), 5.8*radius, 0.8*radius, 0.8*radius, 2*ndiv);
			draw_circle_normal(0.0, 0.8*radius, ndiv, 1, point(tx, ty, 4.0*radius));
			shader.set_specular(0.0, 0.0);
			break;

		case W_PLASMA:
			radius = 0.018;
			shader.set_cur_color(colorRGBA(BLACK, alpha));
			shader.set_specular(0.8, 10.0);
			rot_angle = max(0.0f, 2.0f*(fire_val - 0.7f));
			if (rot_angle != 0.0) fgRotate(rot_angle, -dir.y, dir.x, 0.0);
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
			shader.set_cur_color(colorRGBA(GOLD, alpha));
			draw_circle_normal(0.0075, 0.009, ndiv, 1, 0.15);
			shader.set_cur_color(colorRGBA(BLACK, alpha));
			draw_cylinder(0.15, 0.005, 0.005, 2*ndiv);
			fgPushMatrix();
			fgTranslate(0.0, 0.0, 0.15);
			draw_cylinder(0.07, 0.005, 0.0, 2*ndiv);
			//draw_cylinder(0.18, radius, radius, 2*ndiv);
			{
				colorRGBA color(0.8, 0.6, 1.0, 0.5*alpha);
				if (p_loaded) {set_emissive_only(color, shader);} else {shader.set_cur_color(color);} // emissive when loaded
				fgScale(1.0, 1.0, 0.2);
				for (unsigned i = 0; i < 3; ++i) {draw_sphere_vbo_back_to_front(point(0.0, 0.0, -0.25+0.08*i), 0.01, ndiv, 0);}
				if (p_loaded) {shader.clear_color_e();}
			}
			fgPopMatrix();
			shader.set_specular(0.0, 0.0);
			break;

		case W_M16:
			if ((wmode&1) == 0) { // normal
				radius = 0.0025;
				shader.set_cur_color(colorRGBA(0.04, 0.04, 0.04, alpha));
				shader.set_specular(0.8, 50.0);
				rot_angle = max(0.0, 1.0*fire_val);
				if (rot_angle != 0.0) fgRotate(rot_angle, -dir.y, dir.x, 0.0);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.076), 0.064,     radius,     radius, 2*ndiv, 1);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.000), 0.076, 2.8*radius, 2.0*radius, 2*ndiv, 1);
				draw_cylinder_at(point(0.6*tx, 0.6*ty, 0.136), 0.012, 1.5*radius, 1.5*radius, 2*ndiv, 1);
			}
			else { // shrapnel chaingun
				radius = 0.004;
				float const rdx(1.4*radius*dir.x/rxy), rdy(1.4*radius*dir.y/rxy);
				shader.set_cur_color(colorRGBA(0.12, 0.12, 0.12, alpha));
				shader.set_specular(0.6, 30.0);
				fgTranslate(0.6*tx, 0.6*ty, 0.0);
				fgPushMatrix();
				fgRotate(15.0*rot_counter, 0.0, 0.0, 1.0);
				draw_chaingun_section( rdx,  rdy, radius, ndiv);
				draw_chaingun_section(-rdx, -rdy, radius, ndiv);
				fgRotate(90.0, 0.0, 0.0, 1.0);
				draw_chaingun_section( rdx,  rdy, radius, ndiv);
				draw_chaingun_section(-rdx, -rdy, radius, ndiv);
				shader.set_cur_color(colorRGBA(0.3, 0.3, 0.3, alpha));
				draw_cylinder_at(point(0.0, 0.0, 0.08), 0.004, 2.45*radius, 2.45*radius, 2*ndiv, 1); // outer band
				fgPopMatrix();
			}
			shader.set_specular(0.0, 0.0);
			break;

		case W_SHOTGUN:
			{
				radius = 0.0045;
				float const rdx(radius*dir.x/rxy), rdy(radius*dir.y/rxy);
				shader.set_cur_color(colorRGBA(0.2, 0.2, 0.2, alpha));
				shader.set_specular(0.6, 30.0);
				rot_angle = max(0.0, 8.0*fire_val);
				if (rot_angle != 0.0) fgRotate(rot_angle, -dir.y, dir.x, 0.0);
				fgTranslate(0.6*tx, 0.6*ty, 0.0);
				fgRotate(90.0, 0.0, 0.0, 1.0);
				point const translates[2] = {point(rdx, rdy, 0.0), point(-0.9*rdx, -0.9*rdy, 0.0)};
				
				for (unsigned i = 0; i < 2; ++i) {
					draw_cylinder_at(translates[i], 0.12, radius, radius, 2*ndiv);
					draw_circle_normal(0.0, radius, ndiv, 1, translates[i]+vector3d(0.0, 0.0, 0.09));
				}
				shader.set_specular(0.0, 0.0);
			}
			break;
			
		case W_BBBAT:
			radius = 0.004;
			shader.set_cur_color(colorRGBA(LT_BROWN, alpha));
			select_texture(is_camera ? PLAYER_BBB_TEX : WOOD_TEX); // customize the player's baseball bat
			fgRotate(45.0, -dir.y, dir.x, 0.0);
			fgTranslate(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0); // cancel out texture rotation with camera
			fgRotate(get_bbbat_angle(fire_val), (left_handed ? 0.5 : -0.5), 0.5, 0.0);
			draw_cylinder(16.0*radius, radius, 1.2*radius, ndiv);
			draw_sphere_vbo(point(0.0, 0.0, 16.0*radius), 1.2*radius, ndiv, 1);
			break;

		case W_LASER:
			if (shooter == CAMERA_ID && fired) {
				//beams.push_back(beam3d(0, shooter, (pos0 + dir*0.08), (pos0 + dir*1.08), RED)); // should probably use this instead
				fgPushMatrix();
				set_emissive_only(RED, shader);
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
			shader.set_specular(0.0, 1.0);
			break;

		case W_GASSER:
			radius = 0.14*weapons[W_GASSER].blast_radius;
			shader.set_cur_color(colorRGBA(OLIVE*0.7, alpha));
			shader.set_specular(0.7, 30.0);
			draw_cylinder_at(point(tx, ty, 0.0), 16.0*radius, radius, radius, 2*ndiv);
			draw_circle_normal(0.0, radius, ndiv, 1, point(tx, ty, 8.0*radius));
			shader.set_specular(0.0, 0.0);
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
				set_emissive_only(ORANGE, shader);
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
				set_emissive_only(ORANGE, shader);
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


sphere_t get_weapon_bsphere(int weapon) {

	return sphere_t(); // FIXME: write
}


void draw_weapon_simple(point const &pos, vector3d const &dir, float radius, int cid, int wid, float scale, shader_t &shader) {

	draw_weapon(pos, dir, radius, cid, wid, 0, 0, 0, 1, 0, 2, NO_SOURCE, 0, 1.0, 0.0, 0.0, scale, 0, shader);
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


void draw_weapon_in_hand_real(int shooter, bool draw_pass, shader_t &shader) {

	assert(shooter == CAMERA_ID || shooter < num_smileys);
	assert(sstates != NULL);
	player_state &sstate(sstates[shooter]);
	int const wid(sstate.weapon);
	if (wid == W_UNARMED) return;
	float alpha(1.0);
	float const cradius(object_types[SMILEY].radius);
	bool cull_face(0);

	if (shooter == CAMERA_ID) { // camera/player
		if (camera_view) return;

		if (sstate.powerup == PU_INVISIBILITY) {
			bool const flash(sstate.powerup_time < 4*TICKS_PER_SECOND && ((4*sstate.powerup_time/TICKS_PER_SECOND)&1));
			alpha *= (flash ? 0.8 : 0.4);

			if (wid != W_BLADE) {
				cull_face = 1;
				glEnable(GL_CULL_FACE);
				glCullFace(GL_BACK);
			}
		}
	}
	else { // smiley - out of sync by a frame?
		assert(shooter >= 0 && shooter < num_smileys);
		if (sstate.powerup == PU_INVISIBILITY)         return;
		if (sstate.weapon == W_BALL && game_mode == 2) return; // dodgeball already drawn
	}
	int const cid((shooter == CAMERA_ID) ? camera_coll_id : obj_groups[coll_id[SMILEY]].get_obj(shooter).coll_id);
	vector3d const dir(get_sstate_dir(shooter));
	unsigned const delay(max(1u, weapons[wid].fire_delay));
	float const fire_val((float)sstate.fire_frame/(float)delay);
	point const pos((draw_pass == 0 && wid == W_BLADE) ? sstate.cb_pos : get_sstate_draw_pos(shooter));
	select_texture(WHITE_TEX); // always textured
	draw_weapon(pos, dir, cradius, cid, wid, sstate.wmode, sstate.fire_frame, sstate.plasma_loaded, sstate.p_ammo[wid],
		sstate.rot_counter, delay, shooter, (sstate.cb_hurt > 20), alpha, sstate.dpos, fire_val, 1.0, draw_pass, shader);
	if (shooter == CAMERA_ID) fired = 0;
	if (cull_face) glDisable(GL_CULL_FACE);
}


void draw_weapon_in_hand(int shooter, shader_t &shader) {

	if (!game_mode) return;
	draw_weapon_in_hand_real(shooter, 0, shader);
	scheduled_weapons.insert(shooter); // should not be duplicates, but just in case draw_scheduled_weapons() isn't called
}


void draw_scheduled_weapons() {
	
	if (scheduled_weapons.empty()) return;
	shader_t s;
	setup_smoke_shaders(s, 0.01, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1); // keep_alpha=1; enable smoke?

	for (set<int>::const_iterator i = scheduled_weapons.begin(); i != scheduled_weapons.end(); ++i) {
		draw_weapon_in_hand_real(*i, 1, s);
	}
	s.end_shader();
	scheduled_weapons.clear();
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
	draw_plasma(pos, (pos + pos0), radius, psize, ndiv, 1, 0, shader); // should be sphere mapped texgen?
	select_texture(WHITE_TEX);
	glDisable(GL_CULL_FACE);
	if (psize < 0.9*MAX_PLASMA_SIZE) return;

	// lightning eminating from plasma
	if (is_underwater(spos, 1)) { // under water - suicide
		smiley_collision(shooter, shooter, zero_vector, pos0, PLASMA_LT_DAMAGE, PLASMA_LT_D);
	}
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
		if (p2p_dist(spos, player) < min_dist) min_i = CAMERA_ID; // close enough to hit a the player?
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
		line.draw_lines(); // uses a custom shader
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

	float const scale((world_mode == WMODE_UNIVERSE) ? 0.1 : 1.0); // closer near clip for planets
	float const xy1(0.0006*scale), xy2(0.0002*scale), zval(-0.05*scale);
	float const xy[8] = {-xy1, -xy2, xy1, xy2, 0.0, 0.0, 0.0, 0.0};
	glDisable(GL_DEPTH_TEST);
	point_sprite_drawer psd;

	if (in_zoom) {
		shader_t s;
		s.begin_color_only_shader(color);
		fgScale(2.0, 2.0, 1.0);
		vert_wrap_t verts[8];
		for (unsigned i = 0; i < 8; ++i) {verts[i] = point(xy[i], xy[(i+4)&7], zval);}
		draw_verts(verts, 8, GL_LINES);
		s.end_shader();
	}
	else {
		for (unsigned i = 0; i < 4; ++i) {psd.add_pt(vert_color(point(xy[2*i], xy[(2*i+4)&7], zval), color));}
	}
	psd.add_pt(vert_color(point(0.0, 0.0, zval), color));
	psd.draw(WHITE_TEX, 1.0); // draw with points of size 1 pixel, not blended
	glEnable(GL_DEPTH_TEST);
}


// not sure where this belongs
void draw_teleporters() {

	shader_t s;
	teleporter::shader_setup(s, 4); // RGBA noise
	s.enable();
	s.add_uniform_float("noise_scale", 1.2);
	s.set_cur_color(WHITE);
	enable_blend();

	for (vector<teleporter>::const_iterator i = teleporters.begin(); i != teleporters.end(); ++i) {
		i->draw(s);
	}
	disable_blend();
}


void teleporter::draw(shader_t &s) const {

	float const ACTIVATE_DELAY = 1.0; // in seconds
	float const use_scale(CLIP_TO_01(ACTIVATE_DELAY - (tfticks - last_used_tfticks)/TICKS_PER_SECOND));
	float const draw_radius(get_draw_radius()), light_radius(8.0*draw_radius), use_light_radius(2.0*use_scale*light_radius);
	colorRGBA const c1(blend_color(YELLOW,  RED,    use_scale, 0));
	colorRGBA const c2(blend_color(WHITE,   YELLOW, use_scale, 0));
	colorRGBA const c3(blend_color(LT_BLUE, BLUE,   use_scale, 0));

	if (use_scale > 0.0 && camera_pdu.sphere_visible_test(pos, use_light_radius)) {
		add_dynamic_light(use_light_radius, pos, LT_BLUE);
	}
	if (game_mode && begin_motion && camera_pdu.sphere_visible_test(pos, light_radius)) {
		colorRGBA const lt_color(blend_color(blend_color(c1, c2, fabs(sin(0.05*tfticks)), 0), c3, fabs(cos(0.07*tfticks)), 0));
		add_dynamic_light(light_radius, pos, lt_color);
	}
	if (camera_pdu.sphere_visible_test(pos, draw_radius)) { // draw pos
		s.add_uniform_color("color1i", c1);
		s.add_uniform_color("color1o", c1);
		s.add_uniform_color("color2i", c2);
		s.add_uniform_color("color2o", c2);
		s.add_uniform_color("color3i", c3);
		s.add_uniform_color("color3o", c3);
		s.add_uniform_float("radius", draw_radius);
		s.add_uniform_float("offset", (100.0*pos.x + 0.001*tfticks)); // used as a hash
		s.add_uniform_vector3d("view_dir", (get_camera_pos() - pos).get_norm()); // local object space
		fgPushMatrix();
		translate_to(pos);
		draw_quads();
		fgPopMatrix();

		if (use_scale > 0.9) {
			shader_t cs;
			cs.begin_color_only_shader(colorRGBA(1.0, 1.0, 1.0, 0.5*(use_scale - 0.9)));
			draw_sphere_vbo_back_to_front(pos, draw_radius, N_SPHERE_DIV, 0); // FIXME: use random dissolve texture like exploding starbase
			s.enable();
		}
	}
	if ((display_mode & 0x10) && camera_pdu.sphere_visible_test(dest, 0.25*radius)) { // draw dest (debugging)
		draw_single_colored_sphere(dest, 0.25*radius, N_SPHERE_DIV, BLUE);
		s.enable();
	}
}



