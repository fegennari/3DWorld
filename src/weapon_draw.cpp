// 3D World - Game Mode Weapon Drawing Code
// by Frank Gennari
// 6/24/06

#include "gameplay.h"
#include "physics_objects.h"


vector<int> weap_cobjs;
set<int> scheduled_weapons;

extern bool invalid_ccache, keep_lasers;
extern int game_mode, window_width, window_height, frame_counter, camera_coll_id, display_mode;
extern int num_smileys, left_handed, iticks, do_zoom, camera_view, fired, UNLIMITED_WEAPONS;
extern float fticks, max_proj_rad;
extern obj_type object_types[];
extern obj_group obj_groups[];
extern GLUquadricObj* quadric;
extern vector<spark_t> sparks;
extern vector<laser_beam> lasers;
extern int coll_id[];
extern blood_spot blood_spots[];
extern player_state *sstates;



void laser_beam::draw() const {

	if (shooter == CAMERA_ID) return; // camera (skip for now)
	float const mag(sqrt(intensity));
	colorRGBA c(color);
	c.alpha *= mag;
	draw_line_tquad(pts[0], pts[1], 0.01*mag, 0.01*mag, c, (distant ? ALPHA0 : c), 0);
}


void draw_lasers() {

	if (lasers.empty()) return;
	begin_line_tquad_draw();

	for (unsigned i = 0; i < lasers.size(); ++i) {
		lasers[i].draw();
	}
	if (!keep_lasers) lasers.clear();
	end_line_tquad_draw();
}


void blood_on_camera(unsigned num_spots) {

	num_spots = min(num_spots, NUM_BS/4);
	float const xrand(0.012*((float)window_width/(float)window_height));

	for (unsigned i = 0, n = 0; i < NUM_BS && n < num_spots; ++i) {
		if (blood_spots[i].time <= 0 || (NUM_BS - i <= num_spots - n)) {
			blood_spots[i].pos.assign(rand_uniform(-xrand, xrand), rand_uniform(-0.012, 0.012), -0.02);
			blood_spots[i].size = rand_uniform(3.0, 50.0);
			blood_spots[i].time = int(10 + 3*blood_spots[i].size);
			++n;
		}
	}
}


void show_blood_on_camera() {

	glColor3f(0.7, 0.0, 0.0);
	enable_blend();
	glDisable(GL_LIGHTING);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.5);
	select_texture(BLUR_TEX);
	
	for (unsigned i = 0; i < NUM_BS; ++i) {
		if (blood_spots[i].time > 0) {
			blood_spots[i].time  -= max(1, iticks); // shrink blood spot
			blood_spots[i].size  *= pow(0.99, max(1, iticks));
			blood_spots[i].pos.y -= 0.00001*sqrt(blood_spots[i].size)*fticks;

			if (blood_spots[i].size > 0.1) { // draw it
				glPushMatrix();
				translate_to(blood_spots[i].pos);
				float const size(0.00006*blood_spots[i].size);
				//draw_circle_normal(0.0, 0.3*size, N_SPHERE_DIV, 0);
				//draw_textured_square(size, 0.0, BLUR_TEX);
				draw_tquad(size, size, 0.0, 1);
				glPopMatrix();
			}
			else {
				blood_spots[i].size = 0.0;
				blood_spots[i].time = 0;
			}
		}
	}
	disable_blend();
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_LIGHTING);
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
	bool const DRAW_WEAP_COBJ(0); // for debugging
	int const surfs((wid == W_BLADE || wid == W_M16 || wid == W_SHOTGUN || wid == W_LASER) ? 1 : 0); // no cylinder ends
	cobj_params cp(0.8, BLACK, DRAW_WEAP_COBJ, 1, NULL, weap_cobjs.size(), -1, 1.0, surfs, 0.0, 0.0, 1); // special mode - shadow but no coll
	float rxy, radius;
	vector3d v_trans, dirn(dir*-1); // dir used to be backwards
	point const pos0(get_final_pos(pos, dirn, cradius, 1.0, rxy, v_trans));
	bool const has_xy_comp(dir.x != 0.0 || dir.y != 0.0);

	switch (wid) {
	case W_UNARMED:
		break;

	case W_BALL:
	case W_SBALL:
	case W_LANDMINE:
	case W_GRENADE:
	case W_CGRENADE:
	case W_STAR5:
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
			point pos1(pos0 - dir*cradius), pos2(pos0 - dir*(cradius-dist));
			pos1.z += 0.01;
			pos2.z += 0.01;
			point p1(pos1), p2(p1);
			unsigned const ndiv(max(1U, min(16U, unsigned(fabs(dist)/max_proj_rad + 0.5))));
			vector3d vd((pos2 - pos1)/float(ndiv));

			for (unsigned i = 0; i < ndiv; ++i) { // large and slow, split into smaller cylinders
				p1  = p2;
				p2 += vd;
				weap_cobjs.push_back(add_coll_cylinder(p1, p2, 0.0025, 0.0025, cp));
			}
			vector3d vrot(dir);
			if (has_xy_comp) rotate_vector3d(cross_product(dir, plus_z), PI_TWO, vrot);
			if (fire_val > 0.0) rotate_vector3d(dir, -540.0*fire_val/TO_DEG, vrot);
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
			point pos1(pos + vector3d(pos0, pos)*0.6);
			point pos2(pos1);
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


void update_weapon_cobjs() { // and update cblade and lighting

	clear_weap_cobjs();
	if (!game_mode) return;
	assert(sstates != NULL);
	float const radius(object_types[SMILEY].radius);
	bool const smileys_enabled(obj_groups[coll_id[SMILEY]].enabled);

	for (int i = CAMERA_ID; i < num_smileys; ++i) { // if invisible, don't draw the weapon
		if ((i != CAMERA_ID && !smileys_enabled) || (i == CAMERA_ID && camera_view)) continue;
		player_state &ss(sstates[i]);
		point const pos(get_sstate_draw_pos(i));
		ss.dpos = 0.0;
		
		if (ss.weapon == W_BLADE) {
			ss.cb_pos = pos;
			do_cblade_damage_and_update_pos(ss.cb_pos, i);
		}
		if (ss.powerup == PU_INVISIBILITY) continue; // smiley's/player's shadow will still be there
		float const fire_val(float(ss.fire_frame)/float(max(1U, weapons[ss.weapon].fire_delay)));
		add_weapon_cobj(pos, get_sstate_dir(i), radius, ss.dpos, fire_val, ss.weapon, ss.wmode);
		add_weapon_lights(i);
	}
}


inline void set_color_alpha(colorRGBA color, point const &pos, float alpha, bool shadowed) {

	color.alpha *= alpha;
	set_shadowed_color(color, pos, shadowed);
}


inline void rotate_to_dir(vector3d const &dir, float vadd, float vmult) {

	glRotatef(TO_DEG*vmult*atan2(dir.y, dir.x) + vadd, 0.0, 0.0, 1.0);
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

	rotate_by_vector(dir*-1.0, 0.0); // undo rotation
	rotate_towards_camera(pos);
}


void draw_chaingun_section(float tx, float ty, float radius, int ndiv) {

	glTranslatef(tx, ty, 0.0);
	gluCylinder(quadric, radius, radius, 0.11, 2*ndiv, ndiv);
	draw_circle_normal_at_z(0.08, 0.0, radius, ndiv, 1);
}


void draw_weapon(point const &pos, vector3d dir, float cradius, int cid, int wid, int wmode, int fframe, int p_loaded,
				 int ammo, int rot_counter, unsigned delay, int shooter, int sb_tex, float alpha, float dpos,
				 float fire_val, float scale, bool draw_pass)
{
	bool const just_fired(abs(fframe - (int)delay) <= 1), is_camera(shooter == CAMERA_ID);
	
	if (draw_pass == 1) {
		if (wid == W_PLASMA && p_loaded) {}
		else if ((wid == W_M16 || wid == W_SHOTGUN) && just_fired) {}
		else return; // nothing to draw
	}
	assert(quadric);
	invalid_ccache = 1;
	dir.negate(); // used to be backwards
	float rot_angle, radius, rxy;
	//glDisable(GL_DEPTH_TEST);
	gluQuadricNormals(quadric, GLU_SMOOTH);
	enable_blend();
	set_fill_mode();
	glPushMatrix();
	translate_to(pos);
	uniform_scale(scale);
	rotate_by_vector(dir, 0.0);
	assert(wid < NUM_WEAPONS);
	vector3d v_trans;
	point const pos0(get_final_pos(pos, dir, cradius, scale, rxy, v_trans));
	float const tx(-cradius*dir.x/rxy), ty(-cradius*dir.y/rxy);
	int const light(get_specular_light()); // Note: Doesn't work correctly if more than one light source is enabled
	
	if (draw_pass == 0) { // draw solid objects
		int ndiv(int(get_zoom_scale()*280.0*cradius/(distance_to_camera(pos) + SMALL_NUMBER)));
		ndiv = min(N_SPHERE_DIV, max(3, ndiv));
		unsigned const oid(weapons[wid].obj_id);
		bool do_texture(0);

		if (oid != UNDEF) {
			assert(oid < NUM_TOT_OBJS);
			do_texture = (object_types[oid].tid >= 0);
		}
		int const cobj(weapons[wid].need_weapon ? cid : -1);
		bool const shadowed(!is_visible_to_light_cobj(pos0, light, cradius, cobj, 0));

		switch (wid) {
		case W_UNARMED:
			break;

		case W_BALL:
		case W_SBALL:
		case W_LANDMINE:
		case W_GRENADE: // draw_grenade()?
		case W_CGRENADE:
		case W_STAR5:
			{
				if (do_texture) select_texture((wid == W_BALL) ? select_dodgeball_texture(shooter) : object_types[oid].tid);
				radius = 0.4*object_types[oid].radius;
				if (wid == W_CGRENADE || (wid == W_GRENADE && (wmode & 1))) radius *= 1.2;
				translate_to(v_trans);
				if (do_texture) rotate_to_dir(dir, 90.0, 1.0);  // cancel out texture rotation with camera
				if (do_texture) glRotatef(45.0, 1.0, 0.0, 0.0); // rotate the texture to face the player
				set_color_alpha(object_types[oid].color, pos0, alpha, shadowed);
				if (!shadowed) set_specular(0.8, 40.0);
				draw_sphere_dlist(all_zeros, radius, ndiv, do_texture);
				if (!shadowed) set_specular(0.0, 0.0);
				if (do_texture) glDisable(GL_TEXTURE_2D);
			}
			break;

		case W_BLADE:
			{
				static int lfc(0);
				static float angle(0.0);
				
				if (frame_counter != lfc) {
					lfc    = frame_counter;
					angle += 6.0*fticks;
				}
				radius = 0.032;
				translate_to(v_trans);
				rotate_to_dir(dir, 90.0, 1.0); // cancel out texture rotation with camera
				glPushMatrix();
				glTranslatef(0.0, -0.009, -0.025);
				glScalef(1.0, 1.0, -1.0);
				int const ndiv2(max(4, 2*(is_camera ? N_SPHERE_DIV : ndiv)/3));
				float const len(dpos + radius), dv(4.0*radius/ndiv2);
				assert(dv > 0.0);

				for (float val=0.0; val < len; val += dv) { // draw extendor - large, split into smaller cylinders
					point const cpos(pos0 + dir*(len - val - 0.5*dv)); // center of the segment
					bool const shadowed2(!is_visible_to_light_cobj(cpos, light, cradius, cobj, 0));
					set_color_alpha(colorRGBA(0.49, 0.51, 0.53, 1.0), cpos, alpha, shadowed2);
					draw_fast_cylinder(point(0.0, 0.0, len-val), point(0.0, 0.0, max(0.0f, len-val-dv)), 0.0025, 0.0025, ndiv2, 0, 0);
				}
				glPopMatrix();
				glRotatef(90.0, 1.0, 0.0, 0.0); // into horizontal plane
				glTranslatef(0.0, -0.028, 0.01);
				if (fframe > 0) glRotatef((540.0*fframe)/delay, 0.0, 1.0, 0.0);
				glRotatef(angle, 0.0, 0.0, 1.0);
				set_color_alpha(WHITE, pos0, alpha, shadowed);
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER, 0.95*alpha);
				if (!shadowed) set_specular(0.9, 90.0);
				float dz((ammo > 1) ? -0.025*radius*ammo : 0.0);

				for (int w = 0; w < max(1, ammo); ++w) { // draw a blade for each ammo
					draw_textured_square(radius, dz, (sb_tex ? SAW_B_TEX : SAW_TEX));
					dz += 0.05*radius;
				}
				if (!shadowed) set_specular(0.0, 0.0);
				glDisable(GL_ALPHA_TEST);
				glDisable(GL_TEXTURE_2D);
			}
			break;

		case W_ROCKET:
			radius = 0.95*object_types[ROCKET].radius;
			set_shadowed_color(colorRGBA(0.15, 0.15, 0.15, alpha), pos0, shadowed);
			if (!shadowed) set_specular(0.9, 50.0);
			rot_angle = max(0.0f, 10.0f*(fire_val - 0.7f)); // recoil
			if (rot_angle != 0.0) glRotatef(rot_angle, -dir.y, dir.x, 0.0); // could probably use rotate_into_plus_z()
			glTranslatef(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0);
			gluCylinder(quadric, 0.8*radius, 0.8*radius, 6.8*radius, 2*ndiv, ndiv);
			draw_circle_normal_at_z(5.0*radius, 0.0, 0.8*radius, ndiv, 1);
			// draw the sight
			glTranslatef(0.8*radius, 0.0, 6.5*radius);
			glRotatef(90.0, 0.0, 1.0, 0.0);
			set_shadowed_color((wmode&1) ? BLACK : colorRGBA(0.9, 0.65, 0.05, alpha), pos0, shadowed); // black/gold
			gluCylinder(quadric, 0.15*radius, 0.0, 0.4*radius, ndiv/2, 1);
			if (!shadowed) set_specular(0.0, 0.0);
			break;

		case W_SEEK_D: // similar to rocket
			radius = 0.95*object_types[SEEK_D].radius;
			set_shadowed_color(colorRGBA(0.05, 0.05, 0.05, alpha), pos0, shadowed);
			if (!shadowed) set_specular(0.7, 30.0);
			rot_angle = max(0.0f, 15.0f*(fire_val - 0.8f));
			if (rot_angle != 0.0) glRotatef(rot_angle, -dir.y, dir.x, 0.0);
			glTranslatef(tx, ty, 0.0);
			gluCylinder(quadric, 0.8*radius, 0.8*radius, 5.8*radius, 2*ndiv, ndiv);
			draw_circle_normal_at_z(4.0*radius, 0.0, 0.8*radius, ndiv, 1);
			if (!shadowed) set_specular(0.0, 0.0);
			break;

		case W_PLASMA:
			radius = 0.018;
			set_color_alpha(BLACK, pos0, alpha, shadowed);
			if (!shadowed) set_specular(0.8, 10.0);
			rot_angle = max(0.0f, 2.0f*(fire_val - 0.7f));
			if (rot_angle != 0.0) glRotatef(rot_angle, -dir.y, dir.x, 0.0);
			glTranslatef(0.0, 0.0, 0.15);

			if (shooter != NO_SOURCE) {
				draw_circle_normal(0.0175, 0.018, ndiv, 1);
				rotate_to_dir(dir, 0.0, 1.0);
				glLineWidth(2.0);
				draw_line(point(-0.014, 0.01, 0.0), point(-0.01, -0.001, -0.15));
				glLineWidth(1.0);
				rotate_to_dir(dir, 0.0, -1.0);
			}
			set_color_alpha(RED, pos0, alpha, shadowed);
			draw_circle_normal(0.004, 0.0043, ndiv, 1);
			glTranslatef(0.0, 0.0, -0.15);
			glTranslatef(tx, ty, 0.0);
			set_color_alpha(GOLD, pos0, alpha, shadowed);
			draw_circle_normal_at_z(0.15, 0.0075, 0.009, ndiv, 1);
			set_color_alpha(BLACK, pos0, alpha, shadowed);
			glPushMatrix();
			gluCylinder(quadric, 0.005, 0.005, 0.15, 2*ndiv, ndiv);
			glTranslatef(0.0, 0.0, 0.15);
			gluCylinder(quadric, 0.005, 0.0, 0.07, 2*ndiv, ndiv);
			//gluCylinder(quadric, radius, radius, 0.18, 2*ndiv, ndiv);
			{
				colorRGBA const color(0.8, 0.6, 1.0, 0.5*alpha);
				set_shadowed_color(color, pos0, shadowed);
				if (p_loaded) set_color_e(color);
				glScalef(1.0, 1.0, 0.2);
				draw_sphere_dlist_back_to_front(point(0.0, 0.0, -0.25), 0.01, ndiv, 0, 0);
				draw_sphere_dlist_back_to_front(point(0.0, 0.0, -0.17), 0.01, ndiv, 0, 0);
				draw_sphere_dlist_back_to_front(point(0.0, 0.0, -0.09), 0.01, ndiv, 0, 0);
				if (p_loaded) set_color_e(BLACK);
			}
			glPopMatrix();
			if (!shadowed) set_specular(0.0, 0.0);
			break;

		case W_M16:
			if ((wmode&1) == 0) { // normal
				radius = 0.0025;
				set_shadowed_color(colorRGBA(0.04, 0.04, 0.04, alpha), pos0, shadowed);
				if (!shadowed) set_specular(0.8, 50.0);
				rot_angle = max(0.0, 1.0*fire_val);
				if (rot_angle != 0.0) glRotatef(rot_angle, -dir.y, dir.x, 0.0);
				glTranslatef(0.6*tx, 0.6*ty, 0.0);
				draw_cylinder(point(0.0, 0.0, 0.076), 0.064,     radius,     radius, 2*ndiv, ndiv, 1);
				draw_cylinder(point(0.0, 0.0, 0.000), 0.076, 2.8*radius, 2.0*radius, 2*ndiv, ndiv, 1);
				draw_cylinder(point(0.0, 0.0, 0.136), 0.012, 1.5*radius, 1.5*radius, 2*ndiv, ndiv, 1);
			}
			else { // shrapnel chaingun
				radius = 0.004;
				float const rdx(1.4*radius*dir.x/rxy), rdy(1.4*radius*dir.y/rxy);
				set_shadowed_color(colorRGBA(0.12, 0.12, 0.12, alpha), pos0, shadowed);
				if (!shadowed) set_specular(0.6, 30.0);
				glTranslatef(0.6*tx, 0.6*ty, 0.0);
				glPushMatrix();
				glRotatef(15.0*rot_counter, 0.0, 0.0, 1.0);
				draw_chaingun_section(rdx, rdy, radius, ndiv);
				draw_chaingun_section(-2.0*rdx, -2.0*rdy, radius, ndiv);
				glTranslatef(rdx, rdy, 0.0);
				glRotatef(90.0, 0.0, 0.0, 1.0);
				draw_chaingun_section(rdx, rdy, radius, ndiv);
				draw_chaingun_section(-2.0*rdx, -2.0*rdy, radius, ndiv);
				glTranslatef(rdx, rdy, 0.08);
				set_shadowed_color(colorRGBA(0.3, 0.3, 0.3, alpha), pos0, shadowed);
				draw_cylinder(0.004, 2.45*radius, 2.45*radius, 2*ndiv, 1, 1); // outer band
				glPopMatrix();
			}
			if (!shadowed) set_specular(0.0, 0.0);
			break;

		case W_SHOTGUN:
			{
				radius = 0.0045;
				float const rdx(radius*dir.x/rxy), rdy(radius*dir.y/rxy);
				set_shadowed_color(colorRGBA(0.2, 0.2, 0.2, alpha), pos0, shadowed);
				if (!shadowed) set_specular(0.6, 30.0);
				rot_angle = max(0.0, 8.0*fire_val);
				if (rot_angle != 0.0) glRotatef(rot_angle, -dir.y, dir.x, 0.0);
				glTranslatef(0.6*tx, 0.6*ty, 0.0);
				glRotatef(90.0, 0.0, 0.0, 1.0);
				point const translates[2] = {point(rdx, rdy, 0.0), point(-1.9*rdx, -1.9*rdy, 0.0)};
				
				for (unsigned i = 0; i < 2; ++i) {
					translate_to(translates[i]);
					gluCylinder(quadric, radius, radius, 0.12, 2*ndiv, ndiv);
					draw_circle_normal_at_z(0.09, 0.0, radius, ndiv, 1);
				}
				if (!shadowed) set_specular(0.0, 0.0);
			}
			break;

		case W_BBBAT:
			radius = 0.004;
			set_color_alpha(LT_BROWN, pos0, alpha, shadowed);
			select_texture(WOOD_TEX);
			gluQuadricTexture(quadric, GL_TRUE);
			glRotatef(45.0, -dir.y, dir.x, 0.0);
			glTranslatef(tx, ty, 0.0);
			rotate_to_dir(dir, 0.0, 1.0); // cancel out texture rotation with camera
			glRotatef(get_bbbat_angle(fire_val), (left_handed ? 0.5 : -0.5), 0.5, 0.0);
			gluCylinder(quadric, radius, 1.2*radius, 16.0*radius, ndiv, ndiv);
			draw_sphere_dlist(point(0.0, 0.0, 16.0*radius), 1.2*radius, ndiv, 1, 0);
			gluQuadricTexture(quadric, GL_FALSE);
			glDisable(GL_TEXTURE_2D);
			break;

		case W_LASER:
			set_shadowed_color(colorRGBA(0.45, 0.45, 0.45, alpha), pos0, shadowed);
			if (!shadowed) set_specular(0.8, 50.0);
			glTranslatef(tx, ty, 0.04);
			gluCylinder(quadric, 0.006, 0.0015, 0.16, 2*ndiv, ndiv);
			draw_sphere_dlist(point(0.0, 0.0, 0.0), 0.006, ndiv, 0);
			glTranslatef(0.0, 0.0, -0.04);
			if (!shadowed) set_specular(0.0, 0.0);

			if (shooter == CAMERA_ID && fired) {
				//lasers.push_back(laser_beam(0, shooter, p1, p2)); // should probably use this instead
				set_color(RED);
				set_color_e(RED);
				glTranslatef(-tx, -ty, 0.148);
				glRotatef(30.0, -dir.y, dir.x, 0.0);
				glTranslatef(tx, ty, 0.0);
				gluCylinder(quadric, 0.0007, 0.0007, 0.103, ndiv, 1);
				set_color_e(BLACK);
			}
			break;

		case W_GASSER:
			{
				radius = 0.14*weapons[W_GASSER].blast_radius;
				colorRGBA color(OLIVE*0.7);
				color.alpha = alpha;
				set_shadowed_color(color, pos0, shadowed);
				if (!shadowed) set_specular(0.7, 30.0);
				glTranslatef(tx, ty, 0.0);
				gluCylinder(quadric, radius, radius, 16.0*radius, 2*ndiv, ndiv);
				draw_circle_normal_at_z(8.0*radius, 0.0, radius, ndiv, 1);
				if (!shadowed) set_specular(0.0, 0.0);
			}
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
				glTranslatef(tx, ty, 0.0);
				rotate_to_dir(dir, 0.0, 1.0); // cancel out texture rotation with camera
				draw_plasmaball(pos, shooter);
			}
			break;

		case W_M16:
			if (just_fired) {
				float const size(((wmode&1) == 0) ? 0.02 : 0.0272);
				glTranslatef(0.6*tx, 0.6*ty, 0.0);
				set_color(ORANGE);
				set_color_e(ORANGE);
				glTranslatef(0.0, 0.0, (((wmode&1) == 0) ? 0.15 : 0.12));
				if (!is_camera) rotate_into_camera_dir(pos0, dir); // pos0 is approximate
				draw_textured_square_alpha_test(size, 0.0, BLUR_TEX);
				glDisable(GL_TEXTURE_2D);
				set_color_e(BLACK);
			}
			break;

		case W_SHOTGUN:
			if (just_fired) {
				radius = 0.0042;
				float const rdx(radius*dir.x/rxy), rdy(radius*dir.y/rxy);
				set_color(ORANGE);
				set_color_e(ORANGE);
				glTranslatef(0.6*tx, 0.6*ty, 0.0);
				glRotatef(90.0, 0.0, 0.0, 1.0);
				point const translates[2] = {point(-0.9*rdx, -0.9*rdy, 0.124), point(1.9*rdx, 1.9*rdy, -0.002)};
				
				for (unsigned i = 0; i < 2; ++i) {
					translate_to(translates[i]);
					glPushMatrix();
					if (!is_camera) rotate_into_camera_dir(pos0, dir); // pos0 is approximate
					draw_textured_square_alpha_test(8.0*radius, 0.0, BLUR_TEX); // can't rotate towards camera, already rotated
					glPopMatrix();
				}
				set_color_e(BLACK);
				glDisable(GL_TEXTURE_2D);
			}
			break;
		}
	}
	glPopMatrix();
	disable_blend();
	//glEnable(GL_DEPTH_TEST);
}


void draw_weapon_simple(point const &pos, vector3d const &dir, float radius, int cid, int wid, float scale) {

	draw_weapon(pos, dir, radius, cid, wid, 0, 0, 0, 1, 0, 2, NO_SOURCE, 0, 1.0, 0.0, 0.0, scale, 0);
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


void draw_weapon_in_hand_real(int shooter, bool draw_pass) {

	assert(shooter == CAMERA_ID || shooter < num_smileys);
	assert(sstates != NULL);
	player_state &sstate(sstates[shooter]);
	int const wid(sstate.weapon);
	if (wid == W_UNARMED) return;
	float alpha(1.0);
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
	float const cradius(object_types[SMILEY].radius);
	vector3d const dir(get_sstate_dir(shooter));
	unsigned const delay(max(1u, weapons[wid].fire_delay));
	float const fire_val((float)sstate.fire_frame/(float)delay);
	point const pos((draw_pass == 0 && wid == W_BLADE) ? sstate.cb_pos : get_sstate_draw_pos(shooter));
	draw_weapon(pos, dir, cradius, cid, wid, sstate.wmode, sstate.fire_frame, sstate.plasma_loaded, sstate.p_ammo[wid],
		sstate.rot_counter, delay, shooter, (sstate.cb_hurt > 20), alpha, sstate.dpos, fire_val, 1.0, draw_pass);
	if (shooter == CAMERA_ID) fired = 0;
	if (cull_face) glDisable(GL_CULL_FACE);
}


void draw_weapon_in_hand(int shooter) {

	if (!game_mode) return;
	draw_weapon_in_hand_real(shooter, 0);
	scheduled_weapons.insert(shooter); // should not be duplicates, but just in case draw_scheduled_weapons() isn't called
}


void draw_scheduled_weapons() {

	for (set<int>::const_iterator i = scheduled_weapons.begin(); i != scheduled_weapons.end(); ++i) {
		draw_weapon_in_hand_real(*i, 1);
	}
	scheduled_weapons.clear();
}


void draw_plasmaball(point const &pos0, int shooter) { // and shoot lightning

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
	draw_plasma(pos, radius, psize, ndiv, 0, 1, 0);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_CULL_FACE);
	bool const underwater(is_underwater(spos, 1));
	if (!underwater && (rand()&31) == 0) gen_particles((pos + pos0), 1);
	if (psize < 0.9*MAX_PLASMA_SIZE) return;

	// lightning eminating from plasma
	if (underwater) { // under water - suicide
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

		for (unsigned i = 0; i < 3; ++i) {
			target[i] -= pos[i]*v[i];
		}
		target.normalize();
	}
	set_color(LITN_C);
	set_color_e(LITN_C);

	for (unsigned i = 0; i < unsigned((rand()%6)-2); ++i) { // do we really want to call rand() every time?
		bool const hit(min_i >= CAMERA_ID && rand_float() < 0.6);
		point pos2(pos);

		if (hit) { // targeted to a smiley or the player
			for (i = 0; i < 3; ++i) {
				pos2[i] += rand_uniform(0.9, 1.3)*radius*target[i];
			}
			vadd_rand(pos2, 0.5*radius);
			smiley_collision(min_i, shooter, zero_vector, pos2, PLASMA_LT_DAMAGE, PLASMA_LT_D);
			glLineWidth(2.0);
		}
		else {
			vadd_rand(pos2, 1.7*radius);
		}
		glBegin(GL_LINE_STRIP);
		pos.do_glVertex();
		pos2.do_glVertex();

		for (unsigned j = 0; j < unsigned(rand()%15); ++j) { // do we really want to call rand() every time?
			for (unsigned i = 0; i < 3; ++i) {
				pos2[i] *= rand_uniform(1.01, 1.2);
			}
			pos2.do_glVertex();
		}
		glEnd();
		if (hit) glLineWidth(1.0);
	}
	set_color_e(BLACK);
}


void add_weapon_lights(int shooter) {

	if (sstates == NULL) return;
	assert(shooter >= CAMERA_ID && shooter < num_smileys);
	player_state const &sstate(sstates[shooter]);

	if (sstate.weapon == W_PLASMA && sstate.plasma_loaded) {
		point pos(get_sstate_draw_pos(shooter));
		pos.z += 0.22;
		add_dynamic_light(min(3.5, 45.0*sstate.plasma_size*object_types[PLASMA].radius), pos, get_plasma_color(sstate.plasma_size));
	}
}


void show_crosshair(int do_zoom) {

	float const scale((world_mode == WMODE_UNIVERSE) ? 0.1 : 1.0); // closer near clip for planets
	float const xy1(0.0006*scale), xy2(0.0002*scale), zval(-0.05*scale);
	float const xy[8] = {-xy1, -xy2, xy1, xy2, 0.0, 0.0, 0.0, 0.0};
	glDisable(GL_DEPTH_TEST);
	glColor3f(1.0, 1.0, 1.0);
	enable_blend();
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glDisable(GL_LIGHTING);

	if (do_zoom) {
		glScalef(2.0, 2.0, 1.0);
		glBegin(GL_LINES);

		for (unsigned i = 0; i < 8; ++i) {
			glVertex3f(xy[i], xy[(i+4)&7], zval);
		}
		glEnd();
	}
	else {
		glPointSize(2.0);
		glBegin(GL_POINTS);

		for (unsigned i = 0; i < 8; i += 2) {
			glVertex3f(xy[i], xy[(i+4)&7], zval);
		}
		glEnd();
		glPointSize(1.0);
	}
	draw_point(point(0.0, 0.0, zval));
	disable_blend();
	glDisable(GL_POINT_SMOOTH);
	glDisable(GL_LINE_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}



