// 3D World - OpenGL CS184 Computer Graphics Project - Explosion/blast radius code
// Used in both landscape and universe mode
// by Frank Gennari
// 10/15/05

#include "3DWorld.h"
#include "explosion.h"
#include "ship.h"
#include "ship_util.h"
#include "mesh.h"
#include "gl_ext_arb.h"


bool  const NO_NONETYPE_BRS = 0; // faster but fewer lights
float const EXP_LIGHT_SCALE = 2.0;

int gm_blast(0);
blastr latest_blastr;
vector<blastr> blastrs;
vector<unsigned> available;
vector<explosion> explosions;

extern int iticks, game_mode, display_mode, animate2;
extern GLUquadricObj* quadric;


// duration color1 color2
exp_type_params et_params[NUM_ETYPES] = {
	exp_type_params(1.0, BLACK,  BLACK),    // ETYPE_NONE
	exp_type_params(1.0, YELLOW, RED  ),    // ETYPE_FIRE
	exp_type_params(4.0, WHITE,  BLUE ),    // ETYPE_NUCLEAR
	exp_type_params(0.6, CYAN,   WHITE),    // ETYPE_ENERGY
	exp_type_params(0.7, colorRGBA(0.75, 0.2, 0.4, 1.0), colorRGBA(0.2, 0.0, 0.7, 1.0)), // ETYPE_ATOMIC
	exp_type_params(0.8, GREEN,  LT_GREEN), // ETYPE_PLASMA
	exp_type_params(1.5, colorRGBA(1.0,  0.9, 0.7, 0.4), colorRGBA(0.9, 0.7, 0.3, 0.4)), // ETYPE_EMP
	exp_type_params(2.5, BRONZE_C, DK_RED), // ETYPE_STARB
	exp_type_params(1.7, LT_BLUE, WHITE),   // ETYPE_FUSION
	exp_type_params(1.2, WHITE,   CYAN),    // ETYPE_EBURST
	exp_type_params(0.8, GREEN,   WHITE),   // ETYPE_ESTEAL
	exp_type_params(1.2, WHITE,   WHITE),   // ETYPE_ANIM_FIRE
	exp_type_params(0.8, PURPLE,  LT_BLUE), // ETYPE_SIEGE
};


void explosion::check_pointers() {

	if (source != NULL && !source->is_ok()) source = NULL; // invalid source
	if (parent != NULL && !parent->is_ok()) parent = NULL; // invalid parent
}


void blastr::check_pointers() {

	if (parent != NULL && !parent->is_ok()) parent = NULL; // invalid parent
}


void register_explosion(point const &pos, float radius, float damage, unsigned eflags, int wclass, uobject *src, free_obj const *parent) {

	assert(damage >= 0.0);
	explosions.push_back(explosion(pos, radius, damage, eflags, wclass, src, parent));
}


void apply_explosions() {

	vector<explosion> exps(explosions); // copy so that newly generated explosions are delayed until next frame
	explosions.clear();

	for (unsigned i = 0; i < exps.size(); ++i) {
		apply_explosion(exps[i].pos, exps[i].radius, exps[i].intensity, exps[i].flags, exps[i].wclass, exps[i].source, exps[i].parent);
	}
}


void check_explosion_refs() {

	for (unsigned i = 0; i < explosions.size(); ++i) {
		explosions[i].check_pointers();
	}
	for (unsigned i = 0; i < blastrs.size(); ++i) {
		blastrs[i].check_pointers();
	}
}


void add_blastr(point const &pos, vector3d const &dir, float size, float damage, int time, int src,
				colorRGBA const &color1, colorRGBA const &color2, int type, free_obj const *const parent)
{
	if (NO_NONETYPE_BRS && type == ETYPE_NONE) return;
	assert(size > 0.0 && time > 0);
	blastr br(time, type, src, size, damage, pos, dir, color1, color2, parent);
	br.update();

	if (available.empty()) {
		blastrs.push_back(br);
	}
	else {
		unsigned const ix(available.back());
		available.pop_back();
		assert(ix < blastrs.size());
		blastrs[ix] = br;
	}
	if (type == ETYPE_NUCLEAR && damage > 0.0) { // add flash
		add_blastr(pos, (get_camera_pos() - pos), 1.2*size, 0.0, 2, src, WHITE, WHITE, ETYPE_FUSION, parent);
	}
}


void blastr::update() {

	assert(time > 0 && st_time > 0);
	if (type == ETYPE_ANIM_FIRE) return; // size/color stays constant
	float const cscale(float(time)/float(st_time));
	cur_size = size*(0.3 + 0.5/float(time) + 0.2*(1.0 - cscale));
	blend_color(cur_color, color1, color2, cscale, 1);
	cur_color.alpha *= (0.5 + 0.5*cscale);
}


void blastr::process() const { // land mode

	int const x0(get_xpos(pos.x)), y0(get_ypos(pos.y));
	if (!point_interior_to_mesh(x0, y0)) return;
	colorRGBA light_color(cur_color);
	if (type == ETYPE_ANIM_FIRE) blend_color(light_color, YELLOW, RED, float(time)/float(st_time), 1);
	add_dynamic_light(min(3.5, 4.0*size), pos, light_color);
	if (!animate2) return;
	int ltime(0);

	if (damage > 0.0) {
		float rad(3.0*cur_size/(DX_VAL + DY_VAL));
		int const x1(max(x0 - (int)rad, 1)), x2(min(x0 + (int)rad, MESH_X_SIZE-1));
		int const y1(max(y0 - (int)rad, 1)), y2(min(y0 + (int)rad, MESH_Y_SIZE-1));
		rad *= (DX_VAL + DY_VAL)/SQRT2;
		float const radsq(rad*rad/SQRT2), dscale(2.0E-6*min(2000.0f, damage));

		for (int j = y1; j < y2; ++j) {
			for (int k = x1; k < x2; ++k) {
				point const mpt(get_xval(k), get_yval(j), mesh_height[j][k]);
				float const dist_sq(p2p_dist_sq(pos, mpt));
				if (dist_sq < radsq) surface_damage[j][k] += dscale/(dist_sq + 0.01); // do mesh damage
			}
		}
		//if (time == st_time) // only update grass on the first blast?
		modify_grass_at(pos, 0.5*cur_size, 1, 1, 0, 0); // Note: calling this every time looks better, but is slower
	}
	if (gm_blast == 0 || time < ltime) { // used for object damage
		ltime         = time;
		latest_blastr = *this; // only most recent blast does object damage
	}
	gm_blast = 1;
}


void update_blasts() {

	gm_blast = 0;
	unsigned const nbr(blastrs.size());

	for (unsigned i = 0; i < nbr; ++i) {
		blastr &br(blastrs[i]);
		if (br.time == 0) continue; // blastr not in use

		if (animate2) {
			int const decrement(min(iticks, max(1, (br.st_time >> 1)))); // force it to exist for at least one frame

			if (br.time <= decrement) { // just expired
				br.time = 0;
				available.push_back(i);
				continue;
			}
			br.update();
			br.time -= decrement;
			assert(br.time > 0);
		}
		if (world_mode == WMODE_UNIVERSE) { // universe mode
			float const size(4.0*EXP_LIGHT_SCALE*br.cur_size);
			float const scale(size*size/distance_to_camera_sq(br.pos));

			if (scale > 2.5E-5) {
				if (br.cur_color.alpha > 0.01 /*&& br.type != ETYPE_NONE*/ && univ_sphere_vis_dist(br.pos, 0.2*size)) {
					add_br_light(i, br.pos, size, br.parent);
				}
				if (animate2) add_parts_projs(br.pos, br.cur_size, br.dir, br.cur_color, br.type, br.src, br.parent);
			}
		}
		else if (world_mode == WMODE_GROUND && game_mode && br.damage > 0.0) {
			br.process();
		}
	}
}


void rotate_into_dir(vector3d const &dir, point const &pos) {

	if (dir == zero_vector) {
		rotate_towards_camera(pos);
	}
	else {
		rotate_into_plus_z(dir);
	}
}


void draw_blasts() {

	//RESET_TIME;
	bool const universe(world_mode == WMODE_UNIVERSE);
	glDisable(GL_LIGHTING);
	enable_blend();
	gluQuadricTexture(quadric, GL_TRUE);
	//glEnable(GL_ALPHA_TEST);
	//glAlphaFunc(GL_GREATER, 0.05);

	for (unsigned i = 0; i < blastrs.size(); ++i) {
		blastr const &br(blastrs[i]);
		if (br.time == 0 || br.cur_color.alpha == 0.0)     continue; // expired or transparent
		if (br.type == ETYPE_NONE || br.type == ETYPE_EMP) continue; // no draw
		float const size(br.cur_size);
		point const &pos(br.pos);
		assert(size > 0.0);

		if (universe ? univ_sphere_vis_dist(pos, size) : sphere_in_camera_view(pos, size, 0)) {
			br.cur_color.do_glColor();
			float const timescale(((float)br.time)/(float)br.st_time); // 1.0 to 0.0 => 15 to 0

			if (br.type == ETYPE_ANIM_FIRE) {
				glDepthMask(GL_FALSE);
				select_texture(EXPLOSION_TEX);
				glBegin(GL_QUADS);
				draw_animated_billboard(pos, size, timescale);
				glEnd();
				glDepthMask(GL_TRUE);
				glDisable(GL_TEXTURE_2D);
				continue;
			}
			// use distance_to_camera() for non-universe mode?
			//float const sscale(universe ? 2.2/min(0.02f, distance_to_camera(pos)) : 1.0);
			float const sscale(universe ? 0.4/sqrt(size*distance_to_camera(pos)) : 1.0);
			int const ndiv(max(4, min(N_SPHERE_DIV, int(250.0*size*sscale))));

			switch (br.type) {
			case ETYPE_FIRE:
			case ETYPE_PLASMA:
			case ETYPE_EBURST:
				select_texture(PLASMA_TEX);
				//draw_subdiv_sphere(make_pt_global(pos), size, ndiv, 1, 0); // incorrect bfc due to transforms
				draw_sphere_at_tc(make_pt_global(pos), size, ndiv, 1, 1);
				glDisable(GL_TEXTURE_2D);
				break;

			case ETYPE_ENERGY:
			case ETYPE_ATOMIC:
				{
					select_texture(CLOUD_TEX);
					glEnable(GL_ALPHA_TEST);
					glAlphaFunc(GL_GREATER, 0.4*(1.0 - timescale));
					//glEnable(GL_CULL_FACE);
					glPushMatrix();
					global_translate(pos);
					rotate_about(90.0*timescale, br.dir);
					draw_sphere_dlist(all_zeros, size, ndiv, 1);
					glPopMatrix();
					//glDisable(GL_CULL_FACE);
					glDisable(GL_TEXTURE_2D);
					glDisable(GL_ALPHA_TEST);
				}
				break;

			case ETYPE_FUSION:
			case ETYPE_ESTEAL:
			case ETYPE_STARB:
			case ETYPE_NUCLEAR:
			case ETYPE_SIEGE:
				glDepthMask(GL_FALSE);
				glPushMatrix();
				global_translate(pos);

				if (br.type == ETYPE_STARB) {
					rotate_into_dir((universe ? br.dir : all_zeros), pos);
					select_multitex(BLUR_TEX,  0);
					select_multitex(NOISE_TEX, 1);
					glNormal3f(0.0, 0.0, 1.0);
					glBegin(GL_QUADS);
					draw_one_mult_tex_quad(-2.0*size, -2.0*size, 2.0*size, 2.0*size, 0.0);
					glEnd();
					disable_multitex_a();
				}
				else {
					rotate_into_dir(br.dir, pos);
					draw_textured_quad(2.0*size, 2.0*size, 0.0, BLUR_TEX);
					glDisable(GL_TEXTURE_2D);
				}
				glPopMatrix();
				glDepthMask(GL_TRUE);
				break;

			default:
				assert(0);
			} // switch
		} // visibility test
	} // for i
	//glDisable(GL_ALPHA_TEST);
	gluQuadricTexture(quadric, GL_FALSE);
	disable_blend();
	glEnable(GL_LIGHTING);
	//PRINT_TIME("Draw Blasts");
}


void setup_point_light(point const &pos, colorRGBA const &color, float radius, unsigned gl_light) {

	if (color.alpha == 0.0 || color == BLACK) return; // shouldn't get here

	// set color
	float uambient[4], udiffuse[4];

	for (unsigned d = 0; d < 3; ++d) {
		uambient[d] = 0.2*color[d];
		udiffuse[d] = 1.0*color[d];
	}
	uambient[3] = udiffuse[3] = color[3]; // alpha
	set_colors_and_enable_light(gl_light, uambient, udiffuse);

	// set light attenuation
	assert(radius > 0.0);
	float const atten2(0.1/(EXP_LIGHT_SCALE*radius));
	glLightf(gl_light, GL_CONSTANT_ATTENUATION,  0.5);
	glLightf(gl_light, GL_LINEAR_ATTENUATION,    20.0*atten2);
	glLightf(gl_light, GL_QUADRATIC_ATTENUATION, 5000.0*atten2);
	set_gl_light_pos(gl_light, pos, 1.0); // point light source position
}


bool setup_br_light(unsigned index, point const &pos, unsigned gl_light) {

	assert(index < blastrs.size());
	blastr const &br(blastrs[index]);
	if (br.time == 0) return 0;
	setup_point_light(make_pt_global(br.pos), br.cur_color, br.size, gl_light);
	return 1;
}


bool higher_priority(unsigned first, unsigned second) { // is the priority of <first> greather than <second>?

	assert(first < blastrs.size() && second < blastrs.size());
	return (blastrs[first].priority() > blastrs[second].priority());
}


