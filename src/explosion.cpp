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
#include "shaders.h"
#include "draw_utils.h"


bool  const NO_NONETYPE_BRS = 0; // faster but fewer lights
float const EXP_LIGHT_SCALE = 2.0;

int gm_blast(0);
blastr latest_blastr;
vector<blastr> blastrs;
vector<unsigned> available;
vector<explosion> explosions;

extern int iticks, game_mode, display_mode, animate2;

void calc_lit_uobjects();


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
	if (type == ETYPE_ANIM_FIRE) {br.up_vector = signed_rand_vector_norm();}
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


void blastr::add_as_dynamic_light() const {

	colorRGBA light_color(cur_color);
	if (type == ETYPE_ANIM_FIRE) blend_color(light_color, YELLOW, RED, float(time)/float(st_time), 1);
	add_dynamic_light(min(3.5, 4.0*size), pos, light_color); // Note: 3.5 meant for ground mode, but also acceptable for universe mode (lights are never this large)
}


void blastr::process(int &ltime) const { // land mode

	int const x0(get_xpos(pos.x)), y0(get_ypos(pos.y));
	if (!point_interior_to_mesh(x0, y0)) return;
	add_as_dynamic_light();
	if (!animate2) return;

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
		modify_grass_at(pos, 0.5*cur_size, 1, 1); // crush and burn grass Note: calling this every time looks better, but is slower
	}
	if (gm_blast == 0 || time < ltime) { // used for object damage
		ltime         = time;
		latest_blastr = *this; // only most recent blast does object damage
	}
	gm_blast = 1;
}


void update_blasts() {

	//RESET_TIME;
	gm_blast = 0;
	unsigned const nbr((unsigned)blastrs.size());
	int ltime(0);
	if (world_mode == WMODE_UNIVERSE) {calc_lit_uobjects();}

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
				if (br.cur_color.alpha > 0.01 /*&& br.type != ETYPE_NONE*/ && !is_distant(br.pos, 0.2*size) && univ_sphere_vis(br.pos, size)) {
					add_br_light(i, br.pos, size, br.parent);
				}
				if (animate2) {add_parts_projs(br.pos, br.cur_size, br.dir, br.cur_color, br.type, br.src, br.parent);}
			}
		}
		else if (world_mode == WMODE_GROUND && game_mode && br.damage > 0.0) {
			br.process(ltime);
		}
	} // for i
	//PRINT_TIME("Update Blasts");
}


struct ix_type_pair {
	unsigned ix, type;
	ix_type_pair() {}
	ix_type_pair(unsigned ix_, unsigned type_) : ix(ix_), type(type_) {}
	bool operator<(ix_type_pair const &p) const {return ((type == p.type) ? (ix < p.ix) : (type < p.type));}
};


void draw_blasts() {

	if (blastrs.empty()) return;
	//RESET_TIME;
	bool const universe(world_mode == WMODE_UNIVERSE);
	enable_blend();
	shader_t s;
	s.set_vert_shader("multitex_2");
	s.set_frag_shader("multitex_2");
	s.begin_shader();
	s.add_uniform_int("tex0", 0);
	s.add_uniform_int("tex1", 1);
	//s.add_uniform_float("min_alpha", 0.05);
	select_multitex(WHITE_TEX, 1);
	vector<ix_type_pair> to_draw;

	for (unsigned i = 0; i < blastrs.size(); ++i) {
		blastr const &br(blastrs[i]);
		if (br.time == 0 || br.cur_color.alpha == 0.0)     continue; // expired or transparent
		if (br.type == ETYPE_NONE || br.type == ETYPE_EMP) continue; // no draw
		assert(br.cur_size > 0.0);
		if (!(universe ? univ_sphere_vis_dist(br.pos, br.cur_size) : sphere_in_camera_view(br.pos, br.cur_size, 0))) continue;
		to_draw.push_back(ix_type_pair(i, br.type));
	}
	sort(to_draw.begin(), to_draw.end());
	quad_batch_draw qbd;

	for (vector<ix_type_pair>::const_iterator i = to_draw.begin(); i != to_draw.end(); ++i) {
		blastr const &br(blastrs[i->ix]);
		float const timescale(((float)br.time)/(float)br.st_time);
		bool const begin_type(i == to_draw.begin() || i->type != (i-1)->type);
		bool const end_type  (i+1 == to_draw.end() || i->type != (i+1)->type);

		if (br.type == ETYPE_ANIM_FIRE) {
			if (begin_type) {
				glDepthMask(GL_FALSE);
				select_texture(EXPLOSION_TEX);
			}
			qbd.add_animated_billboard(br.pos, get_camera_pos(), br.up_vector, br.cur_color, br.cur_size, br.cur_size, timescale);
			
			if (end_type) {
				qbd.draw_and_clear();
				glDepthMask(GL_TRUE);
			}
			continue;
		}
		if (begin_type) {set_additive_blend_mode();}

		switch (br.type) {
		case ETYPE_FIRE:
		case ETYPE_PLASMA:
		case ETYPE_EBURST:
			{
				if (begin_type) {select_texture(PLASMA_TEX); glEnable(GL_CULL_FACE);}
				// use distance_to_camera() for non-universe mode?
				//float const sscale(universe ? 2.2/min(0.02f, distance_to_camera(pos)) : 1.0);
				br.cur_color.do_glColor();
				float const sscale(universe ? 0.4/sqrt(br.cur_size*distance_to_camera(br.pos)) : 1.0);
				int const ndiv(max(4, min(N_SPHERE_DIV, int(250.0*br.cur_size*sscale))));
				draw_sphere_vbo(make_pt_global(br.pos), br.cur_size, ndiv, 1);
				if (end_type) {glDisable(GL_CULL_FACE);}
			}
			break;

		case ETYPE_ENERGY:
		case ETYPE_ATOMIC:
			{
				if (begin_type) {
					select_texture(CLOUD_TEX);
					//glEnable(GL_CULL_FACE);
				}
				s.add_uniform_float("min_alpha", 0.4*(1.0 - timescale));
				glPushMatrix();
				global_translate(br.pos);
				rotate_about(90.0*timescale, br.dir);
				br.cur_color.do_glColor();
				float const sscale(universe ? 0.4/sqrt(br.cur_size*distance_to_camera(br.pos)) : 1.0);
				int const ndiv(max(4, min(N_SPHERE_DIV, int(250.0*br.cur_size*sscale))));
				uniform_scale(br.cur_size);
				draw_sphere_vbo_raw(ndiv, 1);
				//draw_sphere_vbo_back_to_front(all_zeros, br.cur_size, ndiv, 1);
				glPopMatrix();
				if (end_type) {/*glEnable(GL_CULL_FACE);*/ s.add_uniform_float("min_alpha", 0.0);}
			}
			break;

		case ETYPE_FUSION:
		case ETYPE_ESTEAL:
		case ETYPE_STARB:
		case ETYPE_NUCLEAR:
		case ETYPE_SIEGE:
			if (begin_type) {
				glDepthMask(GL_FALSE);
				select_multitex(((br.type == ETYPE_FUSION) ? FLARE5_TEX : BLUR_TEX), 0);
				if (br.type == ETYPE_STARB) {select_multitex(NOISE_TEX, 1);}
			}
			if (universe) {
				vector3d const dx(2.0*br.cur_size*cross_product(plus_z, br.dir).get_norm());
				vector3d const dy(2.0*br.cur_size*cross_product(dx,     br.dir).get_norm());
				qbd.add_quad_dirs(make_pt_global(br.pos), dx, dy, br.cur_color);
			}
			else {
				qbd.add_billboard(br.pos, get_camera_pos(), plus_z, br.cur_color, 2.0*br.cur_size, 2.0*br.cur_size);
			}
			if (end_type) {
				qbd.draw_and_clear();
				if (br.type == ETYPE_STARB) {select_multitex(WHITE_TEX, 1);} // set back to white
				glDepthMask(GL_TRUE);
			}
			break;

		default:
			assert(0);
		} // switch
		if (end_type) {set_std_blend_mode();}
	} // for i
	s.end_shader();
	disable_blend();
	//PRINT_TIME("Draw Blasts");
}


void setup_point_light(point const &pos, colorRGBA const &color, float radius, unsigned gl_light) {

	if (color.alpha == 0.0 || color == BLACK) return; // shouldn't get here
	colorRGBA const uambient(color*0.2);
	set_colors_and_enable_light(gl_light, uambient, color);
	assert(radius > 0.0);
	float const atten2(0.1/(EXP_LIGHT_SCALE*radius));
	setup_gl_light_atten(gl_light, 0.5, 20.0*atten2, 5000.0*atten2);
	set_gl_light_pos(gl_light, pos, 1.0); // point light source position
}


bool setup_br_light(unsigned index, point const &pos, unsigned gl_light) {

	assert(gl_light < MAX_GL_LIGHT);
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


