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

vector<blastr> blastrs;
vector<unsigned> available;
vector<explosion> explosions;

extern int iticks, game_mode, display_mode, animate2;
extern float cur_explosion_weight;
extern double tfticks;
extern sphere_t cur_explosion_sphere;

void calc_lit_uobjects();
void draw_one_star(colorRGBA const &colorA, colorRGBA const &colorB, point const &pos, float radius, int ndiv, bool add_halo);


// duration color1 color2
exp_type_params et_params[NUM_ETYPES] = { // for universe mode ship explosions
	exp_type_params(1.0, BLACK,  BLACK),    // ETYPE_NONE
	exp_type_params(1.0, YELLOW, RED  ),    // ETYPE_FIRE
	exp_type_params(4.0, WHITE,  BLUE ),    // ETYPE_NUCLEAR
	exp_type_params(0.6, CYAN,   WHITE),    // ETYPE_ENERGY
	exp_type_params(0.7, colorRGBA(0.75, 0.2, 0.4, 1.0), colorRGBA(0.2, 0.0, 0.7, 1.0)), // ETYPE_ATOMIC
	exp_type_params(0.8, GREEN,  LT_GREEN), // ETYPE_PLASMA
	exp_type_params(1.5, colorRGBA(1.0,  0.9, 0.7, 0.4), colorRGBA(0.9, 0.7, 0.3, 0.4)), // ETYPE_EMP
	exp_type_params(2.5, colorRGBA(1.0,  1.0, 0.6, 1.0), DK_RED), // ETYPE_STARB
	exp_type_params(1.7, LT_BLUE, WHITE),   // ETYPE_FUSION
	exp_type_params(1.2, WHITE,   CYAN),    // ETYPE_EBURST
	exp_type_params(0.8, GREEN,   WHITE),   // ETYPE_ESTEAL
	exp_type_params(1.8, WHITE,   WHITE),   // ETYPE_ANIM_FIRE
	exp_type_params(0.8, PURPLE,  LT_BLUE), // ETYPE_SIEGE
	exp_type_params(1.7, LT_BLUE, WHITE),   // ETYPE_FUSION_ROT
	exp_type_params(2.0, WHITE,   BLACK),   // ETYPE_PART_CLOUD
	exp_type_params(2.0, WHITE,   BLACK),   // ETYPE_PC_ICE
	exp_type_params(3.6, GREEN,   colorRGBA(LT_GREEN, 0.0)), // ETYPE_PBALL
};


void explosion::check_pointers() {
	if (source != NULL && !source->is_ok()) {source = NULL;} // invalid source
	if (parent != NULL && !parent->is_ok()) {parent = NULL;} // invalid parent
}


void blastr::check_pointers() {
	if (parent != NULL && !parent->is_ok()) {parent = NULL;} // invalid parent
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
	for (unsigned i = 0; i < explosions.size(); ++i) {explosions[i].check_pointers();}
	for (unsigned i = 0; i < blastrs.size(); ++i) {blastrs[i].check_pointers();}
}


void add_blastr(point const &pos, vector3d const &dir, float size, float damage, int time, int src, colorRGBA const &color1,
	colorRGBA const &color2, int type, free_obj const *const parent, bool emits_light, float exp_sphere_scale)
{
	if (NO_NONETYPE_BRS && type == ETYPE_NONE) return;
	bool const one_frame_only(time == 0);
	if (one_frame_only) {time = 1;}
	assert(size > 0.0 && time >= 0);
	blastr br(time, type, src, size, damage, pos, dir, color1, color2, parent, one_frame_only, emits_light, exp_sphere_scale);
	br.setup();

	if (available.empty()) {blastrs.push_back(br);}
	else {
		unsigned const ix(available.back());
		available.pop_back();
		assert(ix < blastrs.size());
		blastrs[ix] = br;
	}
	if (type == ETYPE_NUCLEAR && (emits_light || damage > 0.0)) { // add flash
		add_blastr(pos, (get_camera_pos() - pos), 1.2*size, 0.0, 0, src, WHITE, WHITE, ETYPE_FUSION, parent, emits_light, exp_sphere_scale);
	}
}


void blastr::setup() {

	if (type == ETYPE_ANIM_FIRE ) {up_vector = signed_rand_vector_norm();}
	if (type == ETYPE_PART_CLOUD || type == ETYPE_PC_ICE) {cloud_exp.setup(size);}
	update(); // initial update
}

void blastr::update() {

	assert(time > 0 && st_time > 0);
	float const cscale(float(time)/float(st_time));
	if (type != ETYPE_ANIM_FIRE) {cur_size = size*(0.3 + 0.5/float(time) + 0.2*(1.0 - cscale));} // size doesn't increase for anim_fire
	blend_color(cur_color, color1, color2, cscale, 1);
	if (type != ETYPE_ANIM_FIRE) {cur_color.alpha *= (0.5 + 0.5*cscale);}
}

void blastr::add_as_dynamic_light() const {
	add_dynamic_light(min(3.5, 4.0*size), pos, cur_color); // Note: 3.5 meant for ground mode, but also acceptable for universe mode (lights are never this large)
}


void blastr::process() const { // land mode

	float const exp_radius(0.6*size*exp_sphere_scale);

	if (animate2 && exp_radius > 0.0 && camera_pdu.sphere_visible_test(pos, exp_radius)) { // rocket, seekd, cgrenade, raptor
		point const &camera(get_camera_pos());
		float const weight(float(time)/float(st_time));
		float const old_w(cur_explosion_sphere.radius/p2p_dist(camera, cur_explosion_sphere.pos)*cur_explosion_weight);
		float const new_w(exp_radius/p2p_dist(camera, pos)*weight);

		if (cur_explosion_weight == 0.0 || new_w > old_w) { // update if current explosion has higher screen space + intensity influence
			cur_explosion_sphere.pos    = pos;
			cur_explosion_sphere.radius = exp_radius * pow((1.0 - weight), 0.75); // radius expands at 0.75 power over time
			cur_explosion_weight        = weight;
		}
	}
	int const x0(get_xpos(pos.x)), y0(get_ypos(pos.y));
	if (!point_interior_to_mesh(x0, y0)) return;
	if (emits_light) {add_as_dynamic_light();}
	if (!animate2 || damage == 0.0) return;
	float rad(3.0f*cur_size/(DX_VAL + DY_VAL));
	int const irad(ceil(rad));
	int const x1(max(x0 - irad, 1)), x2(min(x0 + irad, MESH_X_SIZE-1));
	int const y1(max(y0 - irad, 1)), y2(min(y0 + irad, MESH_Y_SIZE-1));
	rad *= (DX_VAL + DY_VAL)/SQRT2;
	float const radsq(rad*rad/SQRT2), dscale(2.0E-6*min(2000.0f, damage));

	for (int j = y1; j < y2; ++j) {
		for (int k = x1; k < x2; ++k) {
			point const mpt(get_xval(k), get_yval(j), mesh_height[j][k]);
			float const dist_sq(p2p_dist_sq(pos, mpt));
			if (dist_sq < radsq) {surface_damage[j][k] += dscale/(dist_sq + 0.01);} // do mesh damage
		}
	}
	//if (time == st_time) // only update grass on the first blast?
	modify_grass_at(pos, 0.5*cur_size, 1, 1); // crush and burn grass Note: calling this every time looks better, but is slower
}


bool blastr::next_frame(unsigned i) {

	if (time == 0) return 0; // blastr not in use

	if (animate2) {
		int const decrement(one_frame_only ? (one_frame_seen ? time : 0) : min(iticks, max(1, (st_time >> 1)))); // force it to exist for at least one frame

		if (time <= decrement) { // just expired
			time = 0;
			cloud_exp.clear();
			available.push_back(i);
			return 0;
		}
		update();
		time -= decrement;
		assert(time > 0);
	}
	if (world_mode == WMODE_UNIVERSE) { // universe mode
		float const size(4.0*EXP_LIGHT_SCALE*cur_size), scale(size*size/distance_to_camera_sq(pos));
		if (animate2 && (scale > 2.5E-5 || type == ETYPE_EBURST)) {add_parts_projs(pos, cur_size, dir, cur_color, type, src, parent);}
		if (cur_color.alpha > 0.01 && !is_distant(pos, 0.2*size) && univ_sphere_vis(pos, size)) {add_br_light(i, pos, size, parent);}
	}
	else if (world_mode == WMODE_GROUND && game_mode) {process();}
	one_frame_seen = 1;
	return 1;
}


void update_blasts() {

	//RESET_TIME;
	if (animate2) {
		cur_explosion_sphere = sphere_t(); // reset for next iteration
		cur_explosion_weight = 0.0;
	}
	unsigned const nbr((unsigned)blastrs.size());
	if (world_mode == WMODE_UNIVERSE) {calc_lit_uobjects();}
	for (unsigned i = 0; i < nbr; ++i) {blastrs[i].next_frame(i);}
	//PRINT_TIME("Update Blasts");
}


struct ix_type_pair {
	unsigned ix, type;
	ix_type_pair() : ix(0), type(0) {}
	ix_type_pair(unsigned ix_, unsigned type_) : ix(ix_), type(type_) {}
	bool operator<(ix_type_pair const &p) const {return ((type == p.type) ? (ix < p.ix) : (type < p.type));}
};


/*static*/ void cloud_explosion::draw_setup(vpc_shader_t &s) {
	shader_setup(s, 4); // RGBA noise
	s.enable();
	s.set_cur_color(WHITE);
	//enable_blend();
}
void cloud_explosion::draw(vpc_shader_t &s, point const &pos, float radius, int cloud_exp_type) const {
	//if (!sphere_in_camera_view(pos, radius, 0)) return; // checked by the caller
	switch (cloud_exp_type) {
	case CLOUD_EXP_FIRE: s.set_all_colors(YELLOW, ORANGE,   ORANGE,  RED,   RED,  BLACK); break;
	case CLOUD_EXP_ICE : s.set_all_colors(WHITE,  LT_BLUE,  LT_BLUE, BLUE,  BLUE, BLACK); break;
	default: assert(0);
	}
	s.set_uniform_float(s.rad_loc, radius);
	s.set_uniform_float(s.off_loc, (100.0*sinf(pos.x) + ((world_mode == WMODE_UNIVERSE) ? 0.00005 : 0.0004)*tfticks)); // used as a hash
	s.set_uniform_vector3d(s.vd_loc, (get_camera_pos() - pos).get_norm()); // local object space
	fgPushMatrix();
	global_translate(pos);
	draw_quads();
	fgPopMatrix();
}


bool have_explosions() {

	for (unsigned i = 0; i < blastrs.size(); ++i) {
		blastr const &br(blastrs[i]);
		if (br.time > 0 && br.cur_color.alpha > 0.0 && br.type != ETYPE_NONE) return 1;
	}
	return 0;
}

int get_exp_ndiv(point const &pos, float radius) {
	float const sscale((world_mode == WMODE_UNIVERSE) ? 0.4/max(sqrt(radius*distance_to_camera(pos)), TOLERANCE) : 1.0);
	return max(4, min(N_SPHERE_DIV, int(250.0*radius*sscale)));
}

void draw_billboard_explosion(blastr const &br, quad_batch_draw &qbd, point const &camera, bool aligned) {
	
	if (world_mode == WMODE_UNIVERSE) {
		vector3d const dir(aligned ? (camera - br.pos) : br.dir);
		vector3d const dx(2.0*br.cur_size*cross_product(plus_z, dir).get_norm());
		vector3d const dy(2.0*br.cur_size*cross_product(dx,     dir).get_norm());
		qbd.add_quad_dirs(make_pt_global(br.pos), dx, dy, br.cur_color);
	}
	else {
		qbd.add_billboard(br.pos, camera, plus_z, br.cur_color, 2.0*br.cur_size, 2.0*br.cur_size);
	}
}

void draw_blasts(shader_t &s) {

	if (blastrs.empty()) return;
	//RESET_TIME;
	bool const universe(world_mode == WMODE_UNIVERSE);
	int const min_alpha_loc(s.get_uniform_loc("min_alpha"));
	//s.set_uniform_float(min_alpha_loc, 0.05);
	enable_blend();
	select_texture(WHITE_TEX, 1);
	vector<ix_type_pair> to_draw;
	usw_ray_group exp_rays;

	for (unsigned i = 0; i < blastrs.size(); ++i) {
		blastr const &br(blastrs[i]);
		if (br.time == 0 || br.cur_color.alpha == 0.0)     continue; // expired or transparent
		if (br.type == ETYPE_NONE || br.type == ETYPE_EMP) continue; // no draw
		assert(br.cur_size > 0.0);
		if (!(universe ? univ_sphere_vis_dist(br.pos, br.cur_size) : sphere_in_camera_view(br.pos, br.cur_size, 0))) continue;
		to_draw.push_back(ix_type_pair(i, br.type));
	}
	sort(to_draw.begin(), to_draw.end());
	static quad_batch_draw qbd;
	vpc_shader_t vpc_shader;
	vector<ix_type_pair>::const_iterator range_start(to_draw.begin());
	point const camera(get_camera_pos());

	for (vector<ix_type_pair>::const_iterator i = to_draw.begin(); i != to_draw.end(); ++i) {
		blastr const &br(blastrs[i->ix]);
		float const timescale(((float)br.time)/(float)br.st_time); // goes from 1.0 to 0.0
		bool const begin_type(i == to_draw.begin() || i->type != (i-1)->type);
		bool const end_type  (i+1 == to_draw.end() || i->type != (i+1)->type);
		if (begin_type) {range_start = i;}

		if (br.type == ETYPE_ANIM_FIRE) {
			if (begin_type) {
				glDepthMask(GL_FALSE);
				select_texture(EXPLOSION_TEX);
			}
			float const ts_val(max(0.0, (1.5*timescale - 0.5))); // compact animation into first 67%
			float const alpha(min(1.0, 1.5*timescale)); // drops from 1.0 to 0.0 in second 33%
			qbd.add_animated_billboard(br.pos, camera, br.up_vector, colorRGBA(WHITE, alpha), br.cur_size, br.cur_size, ts_val); // Note: *not* using cur_color
			
			if (end_type) {
				qbd.draw_and_clear();
				glDepthMask(GL_TRUE);
			}
			continue;
		}
		if (begin_type) {set_additive_blend_mode();} // all types below use additive blending

		switch (br.type) {
		case ETYPE_FIRE:
		case ETYPE_PLASMA:
		case ETYPE_EBURST:
			{
				if (begin_type) {select_texture(PLASMA_TEX); glEnable(GL_CULL_FACE); begin_sphere_draw(1);}
				// use distance_to_camera() for non-universe mode?
				//float const sscale(universe ? 2.2/min(0.02f, distance_to_camera(pos)) : 1.0);
				s.set_cur_color(br.cur_color);
				draw_sphere_vbo(make_pt_global(br.pos), br.cur_size, get_exp_ndiv(br.pos, br.cur_size), 1); // cube mapped sphere? too slow?
				
				if (end_type) {
					glDisable(GL_CULL_FACE);
					end_sphere_draw();
					glDepthMask(GL_FALSE);
					select_texture((br.type == ETYPE_FIRE) ? (int)FLARE4_TEX : (int)FLARE2_TEX);

					for (auto j = range_start; j != i+1; ++j) { // go back and add flares
						draw_billboard_explosion(blastrs[j->ix], qbd, camera, 1); // aligned to the camera
					}
					qbd.draw_and_clear();
					glDepthMask(GL_TRUE);
				}
			}
			break;

		case ETYPE_ENERGY:
		case ETYPE_ATOMIC:
			{
				if (begin_type) {select_texture(CLOUD_TEX); /*glEnable(GL_CULL_FACE);*/ begin_sphere_draw(1);}
				s.set_uniform_float(min_alpha_loc, 0.4*(1.0 - timescale));
				fgPushMatrix();
				global_translate(br.pos);
				rotate_about(90.0*timescale, br.dir);
				s.set_cur_color(br.cur_color);
				int const ndiv(get_exp_ndiv(br.pos, br.cur_size));
				uniform_scale(br.cur_size);
				draw_sphere_vbo_raw(ndiv, 1);
				//draw_sphere_vbo_back_to_front(all_zeros, br.cur_size, ndiv, 1);
				fgPopMatrix();
				if (end_type) {/*glEnable(GL_CULL_FACE);*/ end_sphere_draw(); s.set_uniform_float(min_alpha_loc, 0.0);}
			}
			break;

		case ETYPE_FUSION:
		case ETYPE_FUSION_ROT:
		case ETYPE_ESTEAL:
		case ETYPE_STARB:
		case ETYPE_NUCLEAR:
		case ETYPE_SIEGE:
			if (begin_type) {
				glDepthMask(GL_FALSE);
				select_texture(((br.type == ETYPE_FUSION || br.type == ETYPE_FUSION_ROT) ? (int)FLARE5_TEX : ((br.type == ETYPE_STARB) ? (int)STARBURST_TEX : (int)BLUR_TEX)), 0);
			}
			draw_billboard_explosion(br, qbd, camera, (br.type == ETYPE_FUSION_ROT)); // only ETYPE_FUSION_ROT aligns to the camera

			if (end_type) {
				qbd.draw_and_clear();
				glDepthMask(GL_TRUE);
			}
			break;

		case ETYPE_PART_CLOUD:
		case ETYPE_PC_ICE:
			if (begin_type) {cloud_explosion::draw_setup(vpc_shader);}
			vpc_shader.add_uniform_color("color_mult", br.cur_color);
			vpc_shader.add_uniform_float("noise_scale", 0.1/br.size);
			br.cloud_exp.draw(vpc_shader, br.pos, br.cur_size, ((br.type == ETYPE_PC_ICE) ? (int)CLOUD_EXP_ICE : (int)CLOUD_EXP_FIRE)); // Note: no ground mode depth-based attenuation
			
			if (end_type) {
				vpc_shader.add_uniform_color("color_mult", WHITE); // restore default
				vpc_shader.end_shader();
				s.make_current();
			}
			break;

		case ETYPE_PBALL: { // only for universe mode
			rand_gen_t rgen;
			rgen.rseed1 = i->ix;
			unsigned const num(rgen.rand_int(12, 20));

			for (unsigned i = 0; i < num; ++i) { // Note: could use a procedural 1D texture similar to planet rings instead of creating multiple rays
				point const pos2(br.pos + (br.size/SQRT3)*rgen.signed_rand_vector());
				if (pos2 == br.pos) continue; // skip to avoid assertion
				exp_rays.push_back(usw_ray(rgen.rand_uniform(0.01, 0.05)*br.size, rgen.rand_uniform(0.2, 0.4)*br.size, br.pos, pos2, br.cur_color, colorRGBA(br.cur_color, 0.0)));
			}
			if (begin_type) {glEnable(GL_CULL_FACE); select_no_texture();} // texture is procedural
			int const ndiv(get_exp_ndiv(br.pos, br.cur_size));
			colorRGBA color2(br.cur_color, 0.5*br.cur_color.A);
			draw_one_star(br.cur_color, color2, make_pt_global(br.pos), 0.5*br.cur_size, ndiv, 0); // slow, but not used often
			if (end_type) {glDisable(GL_CULL_FACE); exp_rays.draw();}
			break;
		}
		default:
			assert(0);
		} // switch
		if (end_type) {set_std_blend_mode();}
	} // for i
	disable_blend();
	//PRINT_TIME("Draw Blasts");
}

void draw_universe_blasts() {

	if (blastrs.empty()) return;
	shader_t s;
	s.begin_simple_textured_shader();
	draw_blasts(s);
	s.end_shader();
}


void setup_point_light(point const &pos, colorRGBA const &color, float radius, unsigned gl_light, shader_t *shader) {

	if (color.alpha == 0.0 || color == BLACK) return; // shouldn't get here
	colorRGBA const uambient(color*0.2);
	set_colors_and_enable_light(gl_light, uambient, color, shader);
	assert(radius > 0.0);
	float const atten2(0.1f/(EXP_LIGHT_SCALE*radius));
	setup_gl_light_atten(gl_light, 0.5, 20.0*atten2, 5000.0*atten2, shader);
	set_gl_light_pos(gl_light, pos, 1.0, shader); // point light source position
}


bool setup_br_light(unsigned index, point const &pos, unsigned gl_light, shader_t *shader) {

	assert(index < blastrs.size());
	blastr const &br(blastrs[index]);
	if (br.time == 0) return 0;
	setup_point_light(make_pt_global(br.pos), br.cur_color, br.size, gl_light, shader);
	return 1;
}


bool higher_priority(unsigned first, unsigned second) { // is the priority of <first> greather than <second>?

	assert(first < blastrs.size() && second < blastrs.size());
	return (blastrs[first].priority() > blastrs[second].priority());
}


