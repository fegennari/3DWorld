// 3D World - Drawing code for ship/weapon models for universe mode
// by Frank Gennari
// 9/30/05

#include "ship.h"
#include "ship_util.h"
#include "explosion.h"
#include "gl_ext_arb.h"
#include "draw_utils.h"
#include "shaders.h"


bool const TRAIL_FOLLOWS_VEL = 0;
bool const ADD_ENGINE_LIGHTS = 1; // slower and uses lights but looks cool
bool const ADD_CFLASH_LIGHTS = 1; // slower and uses lights but looks cool


extern int display_mode, animate2, frame_counter; // for testing, etc.
extern float fticks;
extern vector<usw_ray> t_wrays;
extern pt_line_drawer_no_lighting_t emissive_pld;


// ******************* USW_RAY, and SHIP_COLL_OBJ classes *******************


void usw_ray::draw(line_tquad_draw_t &drawer) const { // use single sided cylinder with 1D blur rotated towards camera

	// camera view clip?
	//drawer.add_line_tquad(p1, p2, w1, w2, color1, color2, &prev, &next);
	drawer.add_line_as_tris(p1, p2, w1, w2, color1, color2, &prev, &next, 1);
}


void ship_cylinder::draw_cylin(unsigned ndiv, bool textured, float tex_scale_len) const {

	draw_fast_cylinder(p1, p2, r1, r2, ndiv, textured, check_ends, 0, NULL, tex_scale_len);
}


void ship_cube::draw(unsigned ndiv) const { // ndiv is unused

	point pt;
	float sz[3];

	for (unsigned i = 0; i < 3; ++i) {
		pt[i] = 0.5*(d[i][0] + d[i][1]);
		sz[i] =     (d[i][1] - d[i][0]);
	}
	draw_cube(pt, sz[0], sz[1], sz[2], 0);
}


void ship_sphere::draw(unsigned ndiv) const {

	glEnable(GL_CULL_FACE);
	draw_sphere_vbo(pos, radius, ndiv, 0);
	glDisable(GL_CULL_FACE);
}


void ship_torus::draw(unsigned ndiv) const {

	glPushMatrix();
	translate_to(center);
	draw_torus(ri, ro, max(3U, ndiv/2U), ndiv);
	glPopMatrix();
}


void ship_bounded_cylinder::draw(unsigned ndiv) const {
		
	ship_cylinder::draw(ndiv);
	bcube.draw(ndiv);
}


void ship_triangle_list::draw(unsigned ndiv) const { // unused

	vector<vert_norm> verts;
	verts.reserve(3*triangles.size());

	for (vector<triangle>::const_iterator i = triangles.begin(); i != triangles.end(); ++i) {
		vector3d const normal(i->get_normal());
		UNROLL_3X(verts.push_back(vert_norm(i->pts[i_], normal));)
	}
	draw_verts(verts, GL_TRIANGLES);
}


// ******************* GENERAL DRAWING CODE *******************


inline int get_ndiv(int num) {
	int ndiv(max(3, num));
	if (ndiv > 8 && (ndiv&1)) ++ndiv; // an even size is divisible by 2
	return ndiv;
}

void set_ship_texture(int tid) {
	bool const use_shaders((display_mode & 0x08) != 0);
	select_texture(tid, !use_shaders);
}

void end_ship_texture() {
	bool const use_shaders((display_mode & 0x08) != 0);
	end_texture(!use_shaders);
}


void uobj_draw_data::enable_ship_flares(colorRGBA const &color, int tid) {

	set_emissive_color(color);
	glDepthMask(GL_FALSE); // not quite right - prevents flares from interfering with each other but causes later shapes to be drawn on top of the flares
	select_texture(tid);
	set_additive_blend_mode();
}


void uobj_draw_data::disable_ship_flares() {

	qbd.draw_and_clear();
	end_texture();
	glDepthMask(GL_TRUE);
	set_std_blend_mode();
	clear_emissive_color();
}


void setup_colors_draw_flare(point const &pos, point const &xlate, float xsize, float ysize, colorRGBA const &color, int flare_tex) {

	set_emissive_color(color);
	static quad_batch_draw qbd; // probably doesn't need to be static
	qbd.add_xlated_billboard(pos, xlate, get_camera_pos(), up_vector, colorRGBA(0,0,0, color.alpha), xsize, ysize);
	qbd.draw_as_flares_and_clear(flare_tex);
	end_texture();
	clear_emissive_color();
}


void draw_crosshair(upos_point_type const &pos, float dist, colorRGBA const &color) { // for onscreen targeting

	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	enable_blend();
	color.do_glColor();
	glPushMatrix();
	global_translate(pos);
	rotate_towards_camera(pos);
	double const size(0.01*dist);
	vert_wrap_t const verts[3] = {point(0.0, 1.0*size, 0.0), point(-0.7*size, 0.0, 0.0), point( 0.7*size, 0.0, 0.0)}; // z = 0.0
	draw_verts(verts, 3, GL_LINE_LOOP);
	glPopMatrix();
	disable_blend();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
}


void draw_crosshair_from_camera(point const &pos, colorRGBA const &color) {

	float const dist_extend(0.1);
	upos_point_type const &camera(get_player_pos2());
	vector3d_d const v(pos, camera);
	float const vm(v.mag());
	if (vm > TOLERANCE) draw_crosshair((camera + v*(dist_extend/vm)), dist_extend, color);
}


void draw_cobjs(cobj_vector_t const &cobjs, unsigned ndiv) {

	colorRGBA(1.0, 1.0, 1.0, 0.25).do_glColor();

	for (unsigned i = 0; i < cobjs.size(); ++i) {
		assert(cobjs[i]);
		cobjs[i]->draw(ndiv);
	}
}


void us_class::draw_bounding_volume(unsigned ndiv) const {

	if (cobjs.empty()) {
		colorRGBA(1.0, 1.0, 1.0, 0.25).do_glColor();
		bnd_sphere.draw(ndiv); // no special objects (incorrect for dynamic/growing objects)
	}
	else {
		draw_cobjs(cobjs, ndiv);
	}
}


void multipart_ship::draw_bounding_volume(unsigned ndiv) const {

	draw_cobjs(cobjs, ndiv);
}


// ********** uobj_draw_data **********


colorRGBA uobj_draw_data::apply_cloak(colorRGBA const &color) const {

	if (cloakval == 0.0) return color;
	float const cv(CLIP_TO_01(cloakval));
	colorRGBA color2(color);
	color2       *= (1.0 - cv); // scales down both the color components and the alpha value
	color2.alpha *= (1.0 - cv);
	return color2;
}


void uobj_draw_data::draw_bounding_sphere(colorRGBA color) const { // unused

	glEnable(GL_CULL_FACE);
	color.alpha = 0.25;
	color.do_glColor();
	draw_sphere_vbo(all_zeros, crs, ndiv, 0);
	glDisable(GL_CULL_FACE);
}


void uobj_draw_data::setup_exp_scale() const {

	if (t_exp > 0.0) {uniform_scale(1.0 + 0.25*(1.0 - t_exp));} // t_exp drops from 1.0 to 0.0
}


void uobj_draw_data::setup_exp_texture() const {

	if (shader && t_exp > 0.0) { // drops from 1.0 to 0.0
		shader->add_uniform_float("min_alpha", (0.9 + 0.06*(1.0 - t_exp)));
		set_lighted_sides(2);
	}
}


void uobj_draw_data::end_exp_texture() const {

	if (shader && t_exp > 0.0) {
		shader->add_uniform_float("min_alpha", 0.0);
		set_lighted_sides(1);
	}
}


quad_batch_draw uobj_draw_data::qbd;


void uobj_draw_data::draw_engine(colorRGBA const &trail_color, point const &draw_pos, float escale, float ar, vector3d const &stretch_dir) const {

	assert(obj != NULL);
	point viewer(get_player_pos());
	if (ndiv > 3) obj->xform_point(viewer);
	vector3d const orient((viewer - draw_pos).get_norm());
	vector3d up_dir;
	float mod_ar(1.0);

	if (ar == 1.0) {
		up_dir.assign(orient.y, orient.z, orient.x); // swap the xyz values to get an orthogonal vector
	}
	else {
		assert(ar > 1.0 && stretch_dir != all_zeros);
		mod_ar = 1.0 + (ar - 1.0)*cross_product(stretch_dir, orient).mag();
		up_dir = stretch_dir;
	}
	qbd.add_billboard(draw_pos, viewer, up_dir, colorRGBA(0.0, 0.0, 0.0, trail_color.alpha), escale, mod_ar*escale); // color is all emissive

	if (ndiv > 3 && trail_color.alpha != 0.0 && vel.mag_sq() > 1.5E-6) {
		float const dp(dot_product(vel, dir)/vel.mag());
		
		if (dp > 0.0) {
			point epos(draw_pos);
			obj->rotate_point_inv(epos);
			draw_engine_trail(epos, 0.7*escale*sqrt(ar), 0.7, 3.0*dp, trail_color);
		}
	}
}


void uobj_draw_data::draw_engine_trail(point const &offset, float width, float w2s, float len, colorRGBA const &color) const {

	if (!animate2) return;
	if (len < TOLERANCE || (len <= 1.5 && ndiv <= 3) || time < 4) return; // too small/far away to draw
	if (vel.mag_sq() < TOLERANCE*TOLERANCE) return; // not moving
	assert(radius > 0.0 && width > 0.0 && w2s > 0.0 && len > 0.0);
	point const pos2(pos + offset*radius);
	vector3d const delta(len*(TRAIL_FOLLOWS_VEL ? vel : dir*vel.mag())); // 1 tick (not times fticks)
	float const beamwidth(width*radius);
	if (delta.mag_sq() < beamwidth*beamwidth) return; // rarely occurs, but will assertion fail if too small
	t_wrays.push_back(usw_ray(beamwidth, w2s*beamwidth, pos2, (pos2 - delta), color, ALPHA0));
}


void uobj_draw_data::draw_ehousing_pairs(float length, float r1, float r2, float lcone, float dx, float dy, bool texture,
	point const &offset, point const &per_pair_off, unsigned num_pairs) const
{
	assert(length > 0.0 && (r1 > 0.0 || r2 > 0.0));
	unsigned const ndiv2(get_ndiv(ndiv/2));

	for (unsigned p = 0; p < num_pairs; ++p) {
		glPushMatrix();
		translate_to(offset + p*per_pair_off);

		for (unsigned i = 0; i < 2; ++i) { // draw engine housings
			draw_cylin_fast(r1, r2, length, ndiv2, texture);
			if (lcone > 0.0 && r1 > 0.0) draw_cylin_fast(r1, 0.0, lcone, ndiv2, texture);
			glTranslatef(0.0, 0.0, length);
			if (lcone > 0.0 && r2 > 0.0) draw_cylin_fast(r2, 0.0, lcone, ndiv2, texture); // change color?
			if (i == 0) glTranslatef(dx, dy, -length);
		}
		glPopMatrix();
	}
}


void uobj_draw_data::draw_engine_pairs(colorRGBA const &color, unsigned eflags_ix, float escale, float dx, float dy, float dz,
		point const &per_pair_off, unsigned num_pairs, float ar, vector3d const &stretch_dir) const
{
	if (ndiv < 4 || !is_moving() || !first_pass) return; // really should check for thrust, first_pass isn't quite right
	enable_ship_flares(color); // draw engine glow

	for (unsigned p = 0; p < num_pairs; ++p) {
		for (unsigned i = 0; i < 2; ++i) {
			if (!(eflags & (1 << eflags_ix))) {
				draw_engine(color, point((1.0 - 2.0*i)*dx, dy, dz), escale, ar, stretch_dir);
			}
			++eflags_ix;
		}
		dx += per_pair_off.x;
		dy += per_pair_off.y;
		dz += per_pair_off.z;
	}
	disable_ship_flares();
}


bool uobj_draw_data::can_have_engine_lights() const {

	return (ADD_ENGINE_LIGHTS && first_pass && ndiv > 4 && is_moving());
}


void uobj_draw_data::light_engine_pair(colorRGBA const &color, unsigned eflags_off, float escale,
									   float dx, float dy, float dz) const
{
	if (!can_have_engine_lights()) return;
		
	for (unsigned i = 0; i < 2; ++i) {
		if (!(eflags & (1 << (i + eflags_off)))) { // remember that dx and dz are backwards
			setup_point_light(point((1.0 - 2.0*i)*-dx, dy, -dz), color, 4.0*escale*radius, (ENGINE_START_LIGHT + i));
		}
	}
}


void uobj_draw_data::unlight_engine_pair() const {

	if (!can_have_engine_lights()) return;
	for (unsigned i = 0; i < 2; ++i) {clear_colors_and_disable_light(ENGINE_START_LIGHT + i);}
}


void uobj_draw_data::add_light_source(point const &lpos, float lradius, colorRGBA const &color) const {

	if (animate2 && first_pass) add_blastr(lpos, dir, lradius, 0.0, 2, -1, color, color, ETYPE_NONE, obj);
}


void uobj_draw_data::draw_colored_flash(colorRGBA const &color, bool symmetric) const {

	if (!symmetric) glPopMatrix();
	set_emissive_color(color);
	draw_sphere_vbo(all_zeros, 0.25, get_ndiv(ndiv/3), 0); // draw central area that shows up when the draw order is incorrect
	set_emissive_color(color);
	float angle(TWO_PI*time/TICKS_PER_SECOND);
	static quad_batch_draw qbd; // probably doesn't need to be static

	for (unsigned i = 0; i < 2; ++i) {
		qbd.add_xlated_billboard(pos, all_zeros, get_camera_pos(), up_vector, colorRGBA(0,0,0, color.alpha), (2.0 + 1.0*sinf(angle)), (2.0 + 1.0*cosf(angle)));
		angle += PI;
	}
	qbd.draw_as_flares_and_clear(FLARE1_TEX);
	end_texture();
	clear_emissive_color();
	if (!symmetric) glPushMatrix();
	if (ADD_CFLASH_LIGHTS) add_light_source(pos, 4.0*radius, color);
}


void uobj_draw_data::set_cloak_color(colorRGBA const &color) const {

	colorRGBA const cloakc(apply_cloak(color));
	cloakc.do_glColor();
}


void uobj_draw_data::invert_z() const {

	//glScalef(1.0, 1.0, -1.0); // invert z - causes problems with the normals
	glRotatef(180.0, 0.0, 1.0, 0.0); // rotating about y inverts x and z, which isn't quite right, but it doesn't mess up the normals
}


// z is backward, y is up - this is because the default camera in OpenGL looks in -z and the ships use the same coordinate system
void uobj_draw_data::setup_draw_ship() const {

	color_a.do_glColor();
	glPushMatrix();
	invert_z();
}


// ******************* PROJECTILES *******************


void uobj_draw_data::draw_one_triangle(vector3d const &rot_axis, float rot_deg) const {

	vector3d coord_frame[3] = {plus_x, plus_y, -plus_z};
	rotate_vector3d_multi(rot_axis, -rot_deg/TO_DEG, coord_frame, 3); // rotate_about(rot_deg, rot_axis);
	coord_frame[2].do_glNormal(); // using two-sided lighting
	vert_wrap_t const verts[3] = {1.4*coord_frame[1], coord_frame[0], -coord_frame[0]};
	draw_verts(verts, 3, GL_TRIANGLES);
}


void uobj_draw_data::draw_rocket_base(colorRGBA const &cb, colorRGBA const &cn, colorRGBA const &ce,
									  float length, float width, float esize, float tailw) const
{
	cb.do_glColor();
	glTranslatef(0.0, 0.0, 0.5*length);
	draw_cylinder(length, width, width, ndiv, 1, 0, 1);
	
	if (ndiv > 3) {
		cn.do_glColor();
		invert_z(); // draw the correct half
		draw_sphere_vbo(all_zeros, width, min(ndiv, N_SPHERE_DIV), 0, 1);
	}
	glPopMatrix(); // remove the rotations
	glPushMatrix();
	vector3d engine_pos(dir*-(1.5*length + 0.1*esize)); // should already be normalized
	setup_colors_draw_flare(pos, engine_pos, esize, esize, ce, FLARE4_TEX);
	draw_engine_trail(engine_pos, tailw, 0.8, 1.5, ce);
}


void uobj_draw_data::draw_usw_rocket() const {

	draw_rocket_base(LT_GRAY, RED, ORANGE, 2.0, 0.5, 2.0, 0.7);
}


void uobj_draw_data::draw_usw_nukedev() const {

	draw_rocket_base(GRAY, BROWN, BLUE, 1.7, 0.55, 1.8, 0.9);
}


void uobj_draw_data::draw_usw_torpedo() const {

	draw_colored_flash(GREEN, 1);
	draw_engine_trail(all_zeros, 1.0, 1.0, 1.5, GREEN);
}


void uobj_draw_data::draw_spherical_shot(colorRGBA const &color) const {

	if (ndiv <= 3) {
		emissive_pld.add_pt(make_pt_global(pos), color);
		return;
	}
	set_emissive_color(color);
	draw_sphere_vbo(all_zeros, 1.0, min(ndiv, N_SPHERE_DIV/2), 0);
	clear_emissive_color();
}


void uobj_draw_data::draw_usw_energy() const {

	draw_spherical_shot(CYAN);
}


void uobj_draw_data::draw_usw_atomic() const {

	colorRGBA color1(0.5, 0.0, 0.8, 1.0), color2(1.0, 0.4, 0.0, 1.0);
	blend_color(color1, color1, color2, (((float)time)/((float)lifetime)), 0);
	draw_spherical_shot(color1);
}


void uobj_draw_data::draw_usw_emp() const {

	float const alpha(0.5*CLIP_TO_01(1.0f - ((float)time+1)/((float)lifetime+1)));
	glPushMatrix();
	select_texture(SBLUR_TEX);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.01);
	set_emissive_color(colorRGBA(1.0, 0.9, 0.7, alpha));
	glRotatef(90.0, 1.0, 0.0, 0.0);
	draw_sphere_vbo_back_to_front(all_zeros, 1.0, ndiv, 1);
	glDisable(GL_ALPHA_TEST);
	end_texture();
	clear_emissive_color();
	glPopMatrix();
	if (animate2) gen_lightning_from(pos, radius, 0.5, obj);
}


void add_lightning_wray(float width, point const &p1, point const &p2) {

	//t_wrays.push_back(usw_ray(width, width, p1, p2, LITN_C, ALPHA0));
	unsigned const num_segments((rand()&7)+1);
	float const ns_inv(1.0/num_segments);
	point cur(p1), prev(p1);
	vector3d delta((p2 - p1)*ns_inv);
	float const dmag(delta.mag());
	if (dmag < TOLERANCE) return; // shouldn't happen?

	for (unsigned i = 0; i < num_segments; ++i) {
		vadd_rand(delta, 0.25*dmag);
		delta *= dmag/delta.mag(); // normalize
		colorRGBA c[2];
		float w[2];
		
		for (unsigned d = 0; d < 2; ++d) {
			blend_color(c[d], ALPHA0, LITN_C, (i+d)*ns_inv, 1);
			c[d].set_valid_color();
			w[d] = (1.0 - 0.5*(i+d)*ns_inv)*width;
		}
		point const next(cur + delta);
		if (i > 0) {t_wrays.back().next = next;}
		t_wrays.push_back(usw_ray(w[0], w[1], cur, next, c[0], c[1]));
		if (i > 0) {t_wrays.back().prev = prev;}
		prev = cur;
		cur  = next;
		// create recursive forks?
	}
}


void uobj_draw_data::draw_usw_dflare() const {
	
	draw_colored_flash(LT_BLUE, 1);
}


void uobj_draw_data::draw_usw_chaff() const {

	LT_GRAY.do_glColor();
	draw_one_triangle(vector3d(1.0, 1.0, 1.0), PI*time);
}


void uobj_draw_data::draw_usw_rfire() const {

	if (1) {
		if (animate2 && first_pass && (time & 1)) {add_blastr(pos, dir, 1.5*radius, 0.0, 16, -1, WHITE, WHITE, ETYPE_ANIM_FIRE, obj);}
		return;
	}
	float const ctime(CLIP_TO_01(1.0f - ((float)time+1)/((float)lifetime+1)));
	glEnable(GL_CULL_FACE);
	set_emissive_color(colorRGBA(1.0, (0.5 + 0.5*ctime), (0.1 + 0.1*ctime), 1.0));
	select_texture(PLASMA_TEX);
	draw_torus(0.12, 0.72, get_ndiv(ndiv/2), ndiv);
	end_texture();
	glDisable(GL_CULL_FACE);
	clear_emissive_color();
}


void uobj_draw_data::draw_usw_shieldd() const {

	draw_spherical_shot(BLUE);
	draw_engine_trail(all_zeros, 1.5, 1.0, 1.65, BLUE);
}


void uobj_draw_data::draw_usw_thunder() const {

	static float lint(0.5);
	lint += rand_uniform(-0.05, 0.05);
	lint  = CLIP_TO_01(lint);
	float const val(CLIP_TO_01(float(0.5*(rand_uniform(0.47, 0.53) + lint)))); // store last intensity per-thunder?
	float const val2(0.001*(int(100000.0*val) % 999));
	colorRGBA const color(0.1, (0.5 + 0.5*val2), (1.0 - 0.5*val2), 1.0);
	setup_colors_draw_flare(pos, all_zeros, 10.0*val, 10.0*val, color);
	if (ADD_CFLASH_LIGHTS) add_light_source(pos, 4.0*radius, color);
}


void uobj_draw_data::draw_usw_star_int(unsigned ndiv_, point const &lpos, point const &lpos0, float size,
									   float rad, float instability, bool lit) const
{
	colorRGBA const light_color(1.0, 0.95, 0.9);
	set_emissive_color(light_color);
	rotate_about(360.0*rand_float(), signed_rand_vector());
	select_texture(NOISE_TEX);

	mesh2d star_mesh;
	star_mesh.set_size(ndiv_);
	if (instability > 0.0) star_mesh.add_random(instability, -0.5*max(1.0f, instability), instability, 4);
	star_mesh.draw_perturbed_sphere(all_zeros, rad, ndiv_, 1);

	clear_emissive_color();
	glPopMatrix(); // undo transforms
	enable_ship_flares(WHITE);
	draw_engine(ALPHA0, lpos0, 3.0*rad*(1.0 + instability));
	disable_ship_flares();

	if (lit && animate2 && first_pass) {
		add_light_source(lpos, 3.9*size*radius, light_color);
		gen_lightning_from(lpos, 1.5*size*radius, 0.01, obj);
	}
}


void uobj_draw_data::draw_usw_star() const {

	glPushMatrix();
	draw_usw_star_int(ndiv, pos, all_zeros, 1.0, 0.4, 0.5*(float(time+1)/float(lifetime+1)), 1);
}


void uobj_draw_data::draw_usw_seige() const {

	glPopMatrix(); // undo transforms
	glPushMatrix();
	colorRGBA const color(PURPLE);
	setup_colors_draw_flare(pos, all_zeros, 5.0, 5.0, color);

	if (time == 0) {
		add_blastr(pos, dir, 10.0*radius, 0.0, int(0.2*TICKS_PER_SECOND), -1, color, color, ETYPE_NONE, obj);
	}
	else {
		draw_engine_trail(all_zeros, 1.0, 0.1, min(time, 16U), color);
	}
	if (ADD_CFLASH_LIGHTS) add_light_source(pos, 4.0*radius, color);
}


// ******************* SHIPS *******************


void uobj_draw_data::draw_base_fighter(vector3d const &scale) const {
	
	// draw ship sides
	color_a.do_glColor();
	cobj_vector_t const &cobjs(obj->get_cobjs());
	assert(cobjs.size() == 1); // triangles only
	cobjs[0]->draw(0); // ndiv is unused

	// draw engines
	assert(nengines > 1);
	bool const is_scaled(scale != vector3d(1.0, 1.0, 1.0));
	if (is_scaled) glPushMatrix();
	if (is_scaled) glScalef(scale.x, scale.y, scale.z);
	setup_draw_ship();
	float const edy(1.5/(nengines-1));

	if (ndiv > 6) {
		unsigned const ndiv2(get_ndiv(ndiv/2));
		color_b.do_glColor();
		glTranslatef(0.0, -0.75, -2.0);

		for (unsigned i = 0; i < nengines; ++i) {
			draw_cylin_fast(0.3, 0.25, 0.25, ndiv2, 0);
			glTranslatef(0.0, edy, 0.0);
		}
	}
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		point pos2(0.0, -0.75, 2.3);
		enable_ship_flares(YELLOW);
		
		if (ndiv > 4) {
			for (unsigned i = 0; i < nengines; ++i) {
				if (!(eflags & (1 << i))) draw_engine(YELLOW, pos2, 1.0); // this dominates the draw time
				pos2.y += edy;
			}
		}
		else {
			pos2.y += edy;
			if (eflags != 7) draw_engine(YELLOW, pos2, 1.7);
		}
		disable_ship_flares();
	}
	if (is_scaled) glPopMatrix();
}


void uobj_draw_data::draw_xwing() const {

	setup_draw_ship();
	// *** WRITE ***
	glPopMatrix(); // undo invert_z()
}


void uobj_draw_data::draw_us_frigate() const {

	assert(nengines > 1);
	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	glPushMatrix();
	glScalef(1.6, 0.28, 1.2);
	draw_sphere_vbo(all_zeros, 1.0, 3*ndiv/2, textured);
	glPopMatrix();
	color_b.do_glColor();
	glPushMatrix();
	glScalef(0.9, 0.6, 0.8);
	if (specular_en) set_specular(0.5, 80.0);
	draw_sphere_vbo(point(0.0, 0.0, 0.1), 1.0, ndiv, 0); // center
	end_specular();
	glPopMatrix();
	float const edx(1.5/(nengines-1));

	if (ndiv > 3) { // draw engines
		unsigned const ndiv2(get_ndiv(ndiv/2));
		glTranslatef(-0.75, 0.0, -1.5);

		for (unsigned i = 0; i < nengines; ++i) {
			color_b.do_glColor();
			draw_cylin_fast(0.2, 0.22, 1.1, ndiv2, textured);
			
			if (ndiv > 10) {
				(color_b*0.2).do_glColor();
				draw_circle_normal(0.0, 0.202, ndiv2, 1, 0.1);
			}
			glTranslatef(edx, 0.0, 0.0);
		}
	}
	if (textured) end_ship_texture();
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		point pos2(-0.75, 0.0, 1.7);
		enable_ship_flares(YELLOW);
		
		for (unsigned i = 0; i < nengines; ++i) {
			if (!(eflags & (1 << i))) draw_engine(YELLOW, pos2, 1.0);
			pos2.x += edx;
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_us_destroyer() const {

	unsigned const ndiv25(get_ndiv(2*ndiv/5));
	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	glTranslatef(0.0, 0.0, -1.0);
	draw_cylinder(1.7, 0.45, 0.45, ndiv, 1, 1, 0);
	glPushMatrix();
	glScalef(1.0, 1.0, 2.0);
	draw_sphere_vbo(point(0.0, 0.0, 0.85), 0.45, ndiv, textured, 1);
	set_cloak_color(RED);
	
	if (ndiv > 5) { // nose
		if (specular_en) set_specular(0.5, 80.0);
		draw_sphere_vbo(point(0.0, 0.0, 1.13), 0.2, ndiv25, 0, 1);
		end_specular();
	}
	color_b.do_glColor();
	glTranslatef(0.0, 0.0, 0.8);
	draw_cylinder(0.05, 0.7, 0.7, ndiv, 1); // front ring
	glPopMatrix();
	draw_cylindrical_section(point(0.0, 0.0, -0.1), 0.2, 0.85, 1.2, get_ndiv((3*ndiv)/2), textured); // engine ring
	if (textured) end_ship_texture();

	if (ndiv > 4) { // draw engines supports
		draw_cube(all_zeros, 1.8, 0.09, 0.18, 0);
		draw_cube(all_zeros, 0.09, 1.8, 0.18, 0);
	}
	if (ndiv > 6) { // draw parallel beams
		draw_cube(point(0.0, 0.0, 0.8), 1.08,  0.072, 1.62, 0);
		draw_cube(point(0.0, 0.0, 0.8), 0.072, 1.08,  1.62, 0);
	}
	if (ndiv > 8) { // draw engines
		set_cloak_color(GRAY);

		for (unsigned i = 0; i < nengines; ++i) {
			float const theta(TWO_PI*i/((float)nengines));
			glPushMatrix();
			glTranslatef(sinf(theta), cosf(theta), -0.2);
			draw_cylin_fast(0.12, 0.1, 0.14, ndiv25, textured);
			glPopMatrix();
		}
	}
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		float const escale(0.6);
		enable_ship_flares(YELLOW);

		for (unsigned i = 0; i < nengines; ++i) {
			float const theta(TWO_PI*i/((float)nengines));
			if (!(eflags & (1 << i))) draw_engine(YELLOW, point(sinf(theta), cosf(theta), 1.3), escale);
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_us_cruiser(bool heavy) const {

	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	glTranslatef(0.0, 0.0, -0.9);
	draw_cylin_fast(0.38, 0.26, 1.8, ndiv, textured);
	glPushMatrix();
	glScalef(1.0, 1.0, 1.45);
	invert_z(); // invert for half sphere
	draw_sphere_vbo(all_zeros, 0.38, ndiv, textured, 1); // back
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.0, 0.0, 1.8);
	glScalef(1.0, 1.0, 2.12);
	draw_sphere_vbo(all_zeros, 0.26, ndiv, textured, 1); // front
	glPopMatrix();
	float const escale(0.4), sz(2.0*escale);
	point epos[4] = {point(sz, 0, 0.0), point(-2*sz, 0, 0), point(sz, sz, 0), point(0, -2*sz, 0)};
	color_b.do_glColor();
	point pos2(0.0, 0.0, 0.12);
	glPushMatrix();
	glScalef(1.0, 1.0, 4.0);
	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv3(get_ndiv(ndiv/3));

	for (unsigned i = 0; i < nengines; ++i) { // draw engine housings
		pos2 += epos[i];
		draw_sphere_vbo(pos2, 0.16, ndiv2, textured);

		if (ndiv > 8) { // draw engines
			glPushMatrix();
			glTranslatef(pos2.x, pos2.y, (pos2.z - 0.165));
			draw_cylin_fast(0.06, 0.08, 0.03, ndiv3, textured);
			glPopMatrix();
		}
	}
	glPopMatrix();

	if (ndiv > 5) { // draw engine struts
		glPushMatrix();
		point pos0(0.0, 0.0, 0.24);
		glScalef(1.0, 1.0, 2.0);

		for (unsigned i = 0; i < nengines; ++i) {
			pos0   += epos[i];
			point pos2(pos0);
			pos2.x *= 0.2; pos2.y *= 0.2;
			draw_fast_cylinder(pos0, pos2, 0.1, 0.1, ndiv2, textured);
		}
		glPopMatrix();
	}
	if (heavy && ndiv >= 4) { // Note: push/pop not needed since this is the last draw
		float const sz(0.32);
		point const epos[4] = {point(sz, 0, sz), point(-2*sz, 0, 0), point(sz, sz, 0), point(0, -2*sz, 0)};
		point pos2(all_zeros);
		glRotatef(45.0, 0.0, 0.0, 1.0);
		glScalef(1.0, 1.0, 5.0);
		unsigned const ndiv25(get_ndiv((2*ndiv)/5));

		for (unsigned i = 0; i < 4; ++i) { // draw weapon housings
			pos2 += epos[i];
			draw_sphere_vbo(pos2, 0.1, ndiv25, textured);
		}
	}
	if (textured) end_ship_texture();
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		point pos2(0.0, 0.0, 0.46/escale);
		enable_ship_flares(LT_BLUE);
		
		for (unsigned i = 0; i < nengines; ++i) {
			pos2 += epos[i];
			if (!(eflags & (1 << i))) draw_engine(LT_BLUE, pos2, escale);
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_us_bcruiser() const {

	assert(nengines == 4);
	unsigned const ndiv2(get_ndiv(ndiv/2));
	setup_draw_ship();
	float const escale(0.35), dy(0.05), erad(0.1);
	colorRGBA const ecolor(1.0, 0.9, 0.2);
	light_engine_pair(ecolor, 0, escale, 0.4, dy, (1.0 + 0.2*escale)); // average of two engines, eflags isn't quite right
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	// body
	glPushMatrix();
	glTranslatef(0.0, 0.0, -0.8);
	glScalef(0.5, 1.5, 1.5);
	invert_z();
	draw_sphere_vbo(all_zeros, 0.4, ndiv, textured, 1); // rear
	invert_z();
	draw_cylinder(0.6, 0.4, 0.3, ndiv); // main body
	glTranslatef(0.0, 0.0, 0.55);
	glScalef(1.0, 1.0, 3.06);
	draw_sphere_vbo(all_zeros, 0.3, ndiv, textured, 1); // front
	glPopMatrix();

	// engine support wing
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, dy, -0.2);
	glScalef(1.1, 0.18, 1.2);
	draw_sphere_vbo(all_zeros, 0.65, ndiv, textured);
	glPopMatrix();
	set_cloak_color(GRAY);

	if (ndiv > 5) {
		unsigned const ndiv3(get_ndiv(ndiv/3));

		for (unsigned i = 0; i < 2; ++i) { // left, right
			for (unsigned j = 0; j < 2; ++j) { // bottom, top
				point cur(0.08*(2.0*i - 1.0), 0.16*(2.0*j - 1.0), 0.6);
				// draw forward weapons
				if (ndiv > 8) {draw_fast_cylinder(cur, cur+point(0.0, 0.0, 0.7), 0.04, 0.04, ndiv3, textured, ((ndiv > 16) ? 4 : 0));}
				// draw engines
				cur += vector3d(0.32*(2.0*i - 1.0), (dy + (1.5*erad - 0.16)*(2.0*j - 1.0)), -1.4);
				draw_fast_cylinder(cur, cur+point(0.0, 0.0, 0.8), erad, erad, ndiv3, textured, ((ndiv > 12) ? 3 : 0));
				draw_fast_cylinder(cur+point(0.0, 0.0, -0.2), cur, 0.6*erad, erad, ndiv3, textured);
				glPushMatrix();
				translate_to(cur + vector3d(0.0, 0.0, 0.8));
				glScalef(1.0, 1.0, 1.6);
				draw_sphere_vbo(all_zeros, erad, ndiv3, textured, 1);
				glPopMatrix();
			} // for j
		} // for i
	}
	if (textured) end_ship_texture();
	glPopMatrix(); // undo invert_z()

	// draw engine glow (bottom, top)
	draw_engine_pairs(ecolor, 0, escale, 0.4, (dy - 1.5*erad), (1.0 + 0.2*escale), point(0.0, 3.0*erad, 0.0), 2);
	unlight_engine_pair();
}


void uobj_draw_data::draw_us_enforcer() const { // could be better
	
	unsigned const ndiv32(get_ndiv(3*ndiv/2)), ndiv2(get_ndiv(ndiv/2)), nengines(nengines - 1);
	float const epos(0.44);
	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);
	glTranslatef(0.0, 0.0, -1.1);
	draw_cylinder(2.4, 0.65, 0.0, ndiv32, 1); // main body
	color_b.do_glColor();
	if (textured) end_ship_texture();

	if (ndiv > 3) {
		glPushMatrix();
		glTranslatef(0.0, 0.0, -0.09);
		draw_cylin_fast(0.18, 0.16, 0.09, ndiv, 0); // big center engine
		if (ndiv > 5) draw_circle_normal(0.0, 0.17, ndiv32, 0, 0.045);
		glPopMatrix();
	}
	if (ndiv > 6) { // draw engines
		for (unsigned i = 0; i < nengines; ++i) {
			float const theta(TWO_PI*i/((float)nengines));
			glPushMatrix();
			glTranslatef(epos*sinf(theta), epos*cosf(theta), -0.06);
			draw_cylin_fast(0.12, 0.10, 0.06, ndiv2, 0);
			if (ndiv > 8) {draw_circle_normal(0.0, 0.11, ndiv2, 0, 0.03);}
			glPopMatrix();
		}
	}
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		float const escale(0.3);
		colorRGBA const color(1.0, 0.9, 0.2);
		enable_ship_flares(color);
		if (!(eflags & 1)) draw_engine(color, point(0.0, 0.0, 1.26), 2.0*escale); // center engine
		
		for (unsigned i = 0; i < nengines; ++i) {
			float const theta(TWO_PI*i/((float)nengines));
			if (!(eflags & (1 << (i+1)))) draw_engine(color, point(epos*sinf(theta), epos*cosf(theta), 1.21), escale);
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_us_carrier() const {

	unsigned const ndiv2(get_ndiv(ndiv/2));
	assert(nengines == 2);
	setup_draw_ship();
	set_ship_texture(PLASTER_TEX);
	colorRGBA const ecolor(0.2, 0.2, 1.0, 1.0);
	light_engine_pair(ecolor, 0, 0.5, 0.7, 0.0, 1.3);
	draw_cube(point(0.0, 0.0, -0.32), 0.8, 0.38, 2.04, 1, 1, 1.0, 1); // main body

	if (t_exp > 0.0) { // while exploding, the front section breaks off and floats away
		glPushMatrix();

		if (t_exp == 1.0 && animate2 && first_pass) {
			add_blastr((pos + dir*(1.05*radius)), dir, 2.5*radius, 0.0, int(0.8*TICKS_PER_SECOND), -1, WHITE, BLUE, ETYPE_NUCLEAR, obj);
		}
		glTranslatef(0.0, 0.0, 1.1+0.5*(1.0 - t_exp));
		glRotatef(60.0*(1.0 - t_exp), (dir.x+upv.z), (dir.z+upv.x), (dir.y+upv.y));
		glTranslatef(0.0, 0.0, -1.1);
	}
	glPushMatrix();
	glTranslatef(0.0, -0.09, 1.1);
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	draw_cylinder(0.18, 0.2, 0.2, ndiv, 1); // front

	if (ndiv > 6) {
		color_b.do_glColor();
		set_ship_texture(VSTRIPE_TEX);
		glTranslatef(0.0, 0.0, 0.18);
		draw_cylin_fast(0.03, 0.03, 0.1, ndiv, 1);
		if (specular_en) set_specular(0.4, 50.0);
		draw_sphere_vbo(point(0.0, 0.0, 0.1), 0.07, get_ndiv(ndiv/3), 1); // control tower
		end_specular();
	}
	end_ship_texture();
	glPopMatrix();

	color_a.do_glColor();
	vector3d const n[4] = {vector3d(0.9, 0, 0.45), vector3d(-0.9, 0, 0.45), vector3d(0, 0.99, 0.12), vector3d(0, -0.99, 0.12)};
	vert_norm const verts[16] = {
		vert_norm(point( 0.2, -0.09, 1.1), n[0]), vert_norm(point( 0.4, -0.19, 0.7), n[0]), vert_norm(point( 0.4,  0.19, 0.7), n[0]), vert_norm(point( 0.2,  0.09, 1.1), n[0]),
		vert_norm(point(-0.2,  0.09, 1.1), n[1]), vert_norm(point(-0.4,  0.19, 0.7), n[1]), vert_norm(point(-0.4, -0.19, 0.7), n[1]), vert_norm(point(-0.2, -0.09, 1.1), n[1]),
		vert_norm(point( 0.2,  0.09, 1.1), n[2]), vert_norm(point( 0.4,  0.19, 0.7), n[2]), vert_norm(point(-0.4,  0.19, 0.7), n[2]), vert_norm(point(-0.2,  0.09, 1.1), n[2]),
		vert_norm(point(-0.2, -0.09, 1.1), n[3]), vert_norm(point(-0.4, -0.19, 0.7), n[3]), vert_norm(point( 0.4, -0.19, 0.7), n[3]), vert_norm(point( 0.2, -0.09, 1.1), n[3])};
	draw_verts(verts, 16, GL_QUADS);
	if (t_exp > 0.0) glPopMatrix();

	color_b.do_glColor();
	set_ship_texture(SHIP_HULL_TEX);
	draw_ehousing_pairs(1.0, 0.13, 0.14, 0.2, -1.4, 0.0, 1, point(0.7, 0.0, -1.2));
	set_ship_texture(VSTRIPE_TEX);
	glPushMatrix();
	glTranslatef(0.0, 0.19, -0.2);
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	if (ndiv > 3) draw_cylin_fast(0.22, 0.22, 0.1, ndiv2, 1); // weapons turret
	glScalef(1.0, 1.0, 0.5);
	draw_sphere_vbo(point(0.0, 0.0, 0.2), 0.22, ndiv2, 1, 1);
	glPopMatrix();
	end_ship_texture();

	if (ndiv > 9) {
		unsigned const ndiv4(get_ndiv(ndiv/4));
		glPushMatrix();
		glTranslatef(0.0, 0.24, -0.2);

		// draw radar antenna
		glPushMatrix();
		glRotatef(0.8*on_time, 0.0, 1.0, 0.0);
		glRotatef(-30.0, 1.0, 0.0, 0.0);
		draw_cylin_fast(0.008, 0.005, 0.3, ndiv4, 0);
		glScalef(2.0, 1.0, 0.4);
		invert_z();
		if (ndiv > 24) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // wireframe
		draw_sphere_vbo(point(0.0, 0.0, -0.75), 0.06, get_ndiv(ndiv/3), 0, 1);
		if (ndiv > 24) set_fill_mode();
		glPopMatrix();

		// draw energy beam turret
		glTranslatef(0.0, 0.14, 0.0);
		invert_z();
		vector3d turret_dir(tdir);
		turret_dir.y = max(0.0f, turret_dir.y);
		rotate_into_plus_z(turret_dir);
		draw_cylin_fast(0.015, 0.012, 0.4, ndiv4, 0);
		glPopMatrix();
	}
	if (ndiv > 5) draw_cube(point(0.0, 0.0, -0.6), 1.6, 0.09, 0.26, 0); // draw engines support
	glPopMatrix(); // undo invert_z()
	draw_engine_pairs(ecolor, 0, 0.5, 0.7, 0.0, 1.3, all_zeros, 1, 1.5, vector3d(0.0, 1.0, 0.0)); // medium blue, high aspect ratio
	unlight_engine_pair();
}


void uobj_draw_data::draw_armageddon(mesh2d const &surface_mesh) const {

	assert(nengines == 2);
	setup_draw_ship();
	colorRGBA const ecolor(0.7, 0.7, 1.0, 1.0); // lighter blue
	light_engine_pair(ecolor, 0, 0.45, 0.42, -0.62, 0.86);
	
	// main body, with perturb_map and render_map
	(color_a * 0.7).do_glColor(); // helps to reduce light differences when shadowing
	set_ship_texture(NOISE_TEX);
	float const aratio(0.45), xys(aratio*crs);
	glPushMatrix();
	glScalef(xys, xys, crs);
	if (specular_en) set_specular(0.3, 40.0);
	surface_mesh.draw_perturbed_sphere(all_zeros, 1.0, ndiv, 1); // colors change with stencil shadows when textured???
	end_specular();
	glPopMatrix();

	color_b.do_glColor();
	draw_cube(point(0.0, -0.68, 0.0), 0.6, 0.4, 1.2, 1, 1, 1.0, 1); // box
	end_ship_texture();
	unsigned const nbands(5), nspikes(8); // > 1

	if (ndiv > 6) {
		unsigned const ndiv_c(get_ndiv(2*ndiv/3));
		colorRGBA const cgray(apply_cloak(GRAY)), lcgray(apply_cloak(LT_GRAY));
		float const w(0.05), z_start(-0.8*crs), z_step(1.6*crs/(nbands - 1.0));
		glPushMatrix();
		glTranslatef(0.0, 0.0, (z_start - 0.5*w));

		for (unsigned i = 0; i < nbands; ++i) { // draw bands
			float const zpos(i*z_step + z_start), r(fabs(zpos)/crs), radius(xys*sqrt(1.0 - r*r));
			cgray.do_glColor();

			if (ndiv > 9) {
				draw_cylindrical_section(all_zeros, w, 0.9*radius, 1.05*radius, ndiv_c);
			}
			else {
				draw_cylin_fast(1.05*radius, 1.05*radius, w, ndiv_c, 0);
			}
			if (ndiv > 12) { // draw spikes
				lcgray.do_glColor();
				unsigned const ndiv_s(get_ndiv(ndiv/8));
				float const radius2(1.2*(radius + 0.15*xys)), s_width(0.032*(radius + 0.25*xys));

				for (unsigned j = 0; j < nspikes; ++j) {
					float const theta(TWO_PI*j/((float)nspikes)), x(sinf(theta)), y(cosf(theta));
					point const start(radius*x, radius*y, 0.5*w), end(radius2*x, radius2*y, 0.5*w);
					draw_fast_cylinder(start, end, s_width, 0.0, ndiv_s, 0);
				}
			}
			if (i+1 < nbands) glTranslatef(0.0, 0.0, z_step);
		}
		glPopMatrix();
	}
	set_cloak_color(DK_GRAY); // draw engines
	set_ship_texture(SHIP_HULL_TEX);
	draw_ehousing_pairs(0.9, 0.12, 0.12, 0.15, -0.84, 0.0, 1, point(0.42, -0.62, -0.75));
	end_ship_texture();
	glPopMatrix(); // undo invert_z()
	draw_engine_pairs(ecolor, 0, 0.45, 0.42, -0.62, 0.86);
	unlight_engine_pair();
}


void uobj_draw_data::draw_us_shadow() const { // could be improved

	setup_draw_ship();
	set_ship_texture(NOISE_TEX);
	if (cloakval > 0.0) {glEnable(GL_CULL_FACE); glCullFace(GL_BACK);}
	glPushMatrix();
	glScalef(0.7, 0.3, 1.2);
	draw_sphere_vbo(point(0.0, 0.0, 0.2/1.2), 1.0, ndiv, 1);
	glPopMatrix();
	color_b.do_glColor();
	select_texture(CAMOFLAGE_TEX);
	glScalef(0.9, 0.4, 0.8);
	draw_sphere_vbo(point(0.0, 0.0, -0.6/0.8), 1.0, ndiv, 1);
	end_ship_texture();
	if (cloakval > 0.0) glDisable(GL_CULL_FACE);
	glPopMatrix(); // undo invert_z()
}


void uobj_draw_data::draw_defsat() const {

	setup_draw_ship();
	setup_exp_scale();

	// draw main body
	glPushMatrix();
	glTranslatef(0.0, 0.0, -1.3);
	draw_cylinder(1.0, 0.5, 0.5, ndiv, 1); // body
	glTranslatef(0.0, 0.0, 1.0);
	color_b.do_glColor();
	draw_cylin_fast(0.11, 0.0, 1.9, get_ndiv(ndiv/2), 0); // point
	glScalef(1.0, 1.0, 0.6);
	color_a.do_glColor();
	draw_sphere_vbo(point(0.0, 0.0, 1.4), 0.5, ndiv, (t_exp > 0));
	glPopMatrix();

	// draw solar panels
	color_b.do_glColor();
	if (specular_en) set_specular(0.8, 90.0);
	select_texture(PARTB_TEX);
	
	for (unsigned i = 0; i < 2; ++i) {
		draw_cube(point((i ? -1.0 : 1.0), 0.0, 0.0), 1.1, 0.8, 0.1, 1, 1, 4.0);
	}
	end_texture();
	set_cloak_color(GRAY);
	end_specular();
	draw_cube(point(0.0, 0.0, 0.0), 1.1, 0.08, 0.05, 0);
	glPopMatrix(); // undo invert_z()
}


void uobj_draw_data::draw_starbase() const {

	assert(obj);
	int const cyl_ndiv(get_ndiv((2*ndiv)/3)), spoke_ndiv(get_ndiv(ndiv/2));
	cobj_vector_t const &cobjs(obj->get_cobjs());
	assert(cobjs.size() == 8); // should make this more flexible later
	setup_exp_scale();
	set_ship_texture(SPACESHIP1_TEX);
	if (shader && powered) {shader->add_uniform_float("lum_scale", 2.0); shader->add_uniform_float("lum_offset", -1.0);}

	// draw main body (textured)
	WHITE.do_glColor();
	draw_torus(0.2, 1.0, get_ndiv(2*ndiv/3), 4*ndiv/3, 1.0, 8.0); // take from cobjs?

	if (ndiv > 4) { // draw spokes (textured)
		set_cloak_color(color_b);

		for (unsigned i = 2; i < cobjs.size(); ++i) {
			cobjs[i]->draw_cylin(spoke_ndiv, (t_exp > 0.0), 3.0);
		}
	}
	if (shader && powered) {shader->add_uniform_float("lum_scale", 0.0); shader->add_uniform_float("lum_offset", 0.0);}
	end_ship_texture();

	// draw center (team colored)
	color_a.do_glColor();
	cobjs[1]->draw_cylin(cyl_ndiv, (t_exp > 0.0), 2.0);
}


void uobj_draw_data::draw_borg(bool is_cube, bool is_small) const {

	unsigned const ndiv2(is_cube ? max(1, ndiv/2) : get_ndiv((3*ndiv)/2));
	vector3d view_dir(pos, get_player_pos());
	if (is_cube && obj) obj->rotate_point(view_dir);

	if (phase1) {
		select_texture((is_cube && is_small) ? BCUBE_T_TEX : BCUBE_TEX);
		color_b.do_glColor();

		if (is_cube) {
			draw_cube(all_zeros, 1.95, 1.95, 1.95, 1, 0, 1.0, 0, &view_dir);
		}
		else {
			draw_sphere_vbo(all_zeros, 0.97, ndiv2, 1);
		}
	}
	if (phase2) {
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.1);
		select_texture(SMOKE_TEX);
		colorRGBA outer_color(color_a*0.75);
		outer_color.do_glColor();

		if (is_cube) {
			draw_cube(all_zeros, 2.0, 2.0, 2.0, 1, 0, 1.0, 0, &view_dir);
		}
		else {
			draw_sphere_vbo(all_zeros, 1.0, ndiv2, 1);
		}
		glDisable(GL_ALPHA_TEST);
	}
	end_texture();
	
	if (powered && animate2 && final_pass && phase1) {
		add_colored_lights(pos, radius, GREEN, 0.2, (rand() & (is_small ? 3 : 7)), obj); // add eerie green light sources
	}
}


void uobj_draw_data::draw_bshuttle() const {

	setup_draw_ship();
	select_texture(BCUBE_T_TEX);
	draw_cube(point(0.0, 0.0, -0.4), 1.3, 0.6, 1.6, 1);
	end_texture();

	vector3d const n[2] = {vector3d(0.0, 0.8, 0.6), vector3d(0.0, -0.8, 0.6)};
	vert_norm const verts[18] = {
		vert_norm(point( 0.65, 0.0, 1.2), n[0]), vert_norm(point( 0.65,  0.3, 0.4), n[0]), vert_norm(point(-0.65,  0.3, 0.4), n[0]), // T1
		vert_norm(point( 0.65, 0.0, 1.2), n[0]), vert_norm(point(-0.65,  0.3, 0.4), n[0]), vert_norm(point(-0.65,  0.0, 1.2), n[0]), // T2
		vert_norm(point(-0.65, 0.0, 1.2), n[1]), vert_norm(point(-0.65, -0.3, 0.4), n[1]), vert_norm(point( 0.65, -0.3, 0.4), n[1]), // B1
		vert_norm(point(-0.65, 0.0, 1.2), n[1]), vert_norm(point( 0.65, -0.3, 0.4), n[1]), vert_norm(point( 0.65,  0.0, 1.2), n[1]), // B2
		vert_norm(point(-0.65, 0.3, 0.4), -plus_x), vert_norm(point(-0.65, -0.3, 0.4), -plus_x), vert_norm(point(-0.65, 0.0, 1.2), -plus_x),  // L
		vert_norm(point( 0.65, 0.0, 1.2),  plus_x), vert_norm(point( 0.65, -0.3, 0.4),  plus_x), vert_norm(point( 0.65, 0.3, 0.4),  plus_x)}; // R
	draw_verts(verts, 18, GL_TRIANGLES);

	color_b.do_glColor();
	select_texture(BCUBE_T_TEX);
	draw_ehousing_pairs(1.0, 0.18, 0.18, 0.22, 1.6, 0.0, 1, point(-0.8, 0.0, -0.8)); // length r1 r2 lcone dx dy
	end_texture();
	glPopMatrix(); // undo invert_z()
	draw_engine_pairs(LT_BLUE, 0, 0.8, 0.8, 0.0, 1.0); // escale dx dy dz
}


void uobj_draw_data::draw_tractor() const { // could be better

	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);
	setup_draw_ship();
	draw_cube(point(0.0, 0.0, -0.2), 1.6, 1.2, 2.0, textured);
	color_b.do_glColor();
	draw_cube(point(0.0, 0.0, 1.0), 0.4, 0.2, 0.4, textured);
	draw_ehousing_pairs(1.0, 0.25, 0.25, 0.3, 2.1, 0.0, 1, point(-1.05, 0.25, -0.8), point(0.0, -0.5, 0.0), 2); // length r1 r2 lcone dx dy
	if (textured) end_ship_texture();
	set_cloak_color(colorRGBA(1.0, 1.0, 1.0, 0.5));
	draw_cube(point(0.0, 0.0, 1.1), 0.8, 0.4, 0.65, 0);
	glPopMatrix(); // undo invert_z()
	draw_engine_pairs(WHITE, 0, 0.9, 1.0, 0.25, 1.0, point(0.0, -0.5, 0.0), 2); // escale dx dy dz
}


void uobj_draw_data::draw_gunship() const {

	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv3(get_ndiv(ndiv/3)), ndiv4(get_ndiv(ndiv/4));
	setup_draw_ship();
	uniform_scale(1.125);

	// front
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);
	glPushMatrix();
	glScalef(1.0, 1.0, 2.0);
	draw_sphere_vbo(point(0.0, 0.0, 0.325), 0.25, 2*ndiv/3, textured);
	glPopMatrix();

	// body
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.0, -0.95);
	glScalef(0.25, 1.0, 1.0);
	draw_cylinder(1.6, 0.45, 0.25, ndiv2); // main body
	glScalef(4.0, 0.25, 1.0);
	draw_cylinder(1.6, 0.45, 0.25, ndiv2); // main body
	glTranslatef(0.0, 0.0, -0.2);
	draw_cylinder(0.2, 0.00, 0.45, ndiv2); // rear
	glScalef(0.25, 4.0, 1.0);
	draw_cylinder(0.2, 0.00, 0.45, ndiv2); // rear
	glPopMatrix();
	if (textured) end_ship_texture();

	// rings
	GOLD.do_glColor();
	glPushMatrix();
	glScalef(3.0, 1.0, 1.0);
	glTranslatef(0.0, 0.0, -0.7);

	for (unsigned i = 0; i < 5; ++i) {
		draw_torus(0.02, (0.3 - 0.05*abs(int(i) - 2)), ndiv3, ndiv);
		glTranslatef(0.0, 0.0, 0.3);
	}
	glPopMatrix();

	// engines
	float const dxy[2][4] = {{-1.0, 1.0, 0.0, 0.0}, {0.0, 0.0, -1.0, 1.0}};

	if (ndiv > 4) { // engine housings
		GRAY.do_glColor();

		for (unsigned i = 0; i < 4; ++i) {
			glPushMatrix();
			glTranslatef(0.35*dxy[0][i], 0.35*dxy[1][i], -1.1);
			draw_cylinder(0.25, 0.05, 0.05, ndiv2);
			if (ndiv > 12) {draw_circle_normal(0.0, 0.05, ndiv2, 0, 0.05);}
			glPopMatrix();
		}
	}
	glPopMatrix(); // undo invert_z()

	if (is_moving()) { // draw engine glow
		enable_ship_flares(GOLD);

		for (unsigned i = 0; i < 4; ++i) {
			point const epos(0.35*dxy[0][i], 0.35*dxy[1][i], 1.15);
			if (!(eflags & (1 << i))) draw_engine(GOLD, epos*1.125, 0.5);
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_nightmare() const {

	unsigned const ndiv23(get_ndiv(2*ndiv/3));
	setup_draw_ship();
	uniform_scale(0.8);
	draw_sphere_vbo(point(0.0, 0.0, -0.1), 1.0, 3*ndiv/2, 0);
	color_b.do_glColor();
	draw_sphere_vbo(point(0.0, 0.0,  0.1), 1.0, 3*ndiv/2, 0);

	if (ndiv > 3) {
		// draw points
		for (unsigned i = 0; i < 4; ++i) {
			float const x(1.0 - 2.0*(i&1)), y(1.0 - 2.0*(i>>1));
			point const pt1(0.5*x, 0.5*y, 0.5), pt2(x, y, 1.0);
			draw_fast_cylinder(pt1, pt2, 0.45, 0.0, 4, 0); // 4 sided
		}
		draw_fast_cylinder(point(0.0, 0.0, 0.8), point(0.0, 0.0, 1.6), 0.45, 0.0, ndiv, 0);
	
		GRAY.do_glColor();
		glTranslatef(0.0, 0.0, -1.2);
		draw_cylin_fast(0.4, 0.35, 0.2, ndiv23, 0); // engine
		draw_circle_normal(0.0, 0.4, ndiv23, 0);
	}
	glPopMatrix(); // undo invert_z()

	if (is_moving() && !(eflags & 1)) { // draw engine glow
		enable_ship_flares(RED);
		draw_engine(RED, point(0.0, 0.0, 1.3), 1.5);
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_dwcarrier() const {

	unsigned const ndiv35(get_ndiv(3*ndiv/5)), ndiv2(get_ndiv(ndiv/2)), ndiv4(get_ndiv(ndiv/4));
	setup_draw_ship();
	if (powered && first_pass) setup_point_light(point(0.0, 0.4, 0.8), RED, 2.0*radius, ENGINE_DEF_LIGHT);

	if (phase1) {
		if (ndiv > 3) { // command center
			glPushMatrix();
			glScalef(1.0, 1.2, 2.8);
			draw_sphere_vbo(point(0.0, 0.22, -0.2), 0.08, ndiv2, 0);
			glPopMatrix();
		}
		color_b.do_glColor();
		glPushMatrix();
		glTranslatef(0.0, 0.0, -1.1);
		
		for (unsigned i = 0; i < 4; ++i) { // draw body
			if (ndiv < 16 && i == 1) i = 3;
			glPushMatrix();
			glScalef((1.0 - 0.18*i), (0.16 + 0.1*i), 1.0);

			if (ndiv < 10) { // approximate low-detail version with only one cylinder and no ends
				draw_cylin_fast(0.72, 0.15, 2.4, ndiv35, 0);
			}
			else {
				glPushMatrix();

				for (unsigned j = 0; j < 4; ++j) {
					float const length(0.93 - 0.2*j);
					draw_cylinder(length, (0.72 - 0.2*j), (0.61 - 0.2*j), ndiv35, (ndiv >= 8), 0, 1);
					if (j < 3) glTranslatef(0.0, 0.0, length);
				}
				glPopMatrix();
			}
			glScalef(1.0, 1.0, 0.36);
			invert_z();
			draw_sphere_vbo(all_zeros, 0.72, ndiv35, 0, 1); // rear
			glPopMatrix();
			glTranslatef(0.0, 0.0, 0.01);
		}
		glPopMatrix();

		if (ndiv > 4) { // engine housings
			glPushMatrix();
			glScalef(0.32, 1.0, 1.0);
			glTranslatef(1.0, 0.0, -1.35);

			for (unsigned i = 0; i < 2; ++i) {
				draw_cylinder(1.0, 0.38, 0.1, ndiv2, (ndiv >= 8), 1, 0);
				if (i == 0) glTranslatef(-2.0, 0.0, 0.0);
			}
			glPopMatrix();
		}
		if (ndiv > 5) { // "eyes"
			RED.do_glColor();
			glPushMatrix();
			glScalef(1.0, 1.0, 2.0);

			for (unsigned i = 0; i < 2; ++i) {
				draw_sphere_vbo(point(0.06*(1.0 - 2.0*i), 0.1, 0.4), 0.025, ndiv4, 0);
			}
			glPopMatrix();
		}
	} // end phase1
	if (phase2 && ndiv > 6) { // razorback
		if (ndiv > 9) { // draw "holes" in backwards depth order
			ALPHA0.do_glColor();
			glPushMatrix();
			glRotatef(90.0, 0.0, 1.0, 0.0);
			glTranslatef(0.84, 0.44, -0.01);

			for (unsigned i = 0; i < 5; ++i) {
				draw_cylinder(0.02, 0.05, 0.05, ndiv4, 1);
				if (i < 4) glTranslatef(-0.3, -0.04, 0.0);
			}
			glPopMatrix();
		}
		color_b.do_glColor();
		glRotatef(10.0, 1.0, 0.0, 0.0); // Note: push/pop not needed since this is the last draw
		glScalef(0.01, 0.2, 1.0);
		draw_sphere_vbo(point(0.0, 0.7, -0.1), 1.0, ndiv4, 0);
	}
	glPopMatrix(); // undo invert_z()

	if (phase2) { // draw engine glow
		draw_engine_pairs(BLUE, 0, 0.4, 0.32, 0.3, 1.45, point(0.0, -0.3, 0.0), 3);
		//draw_engine_pairs(BLUE, 0, 0.4, 0.32, 0.0, 1.45, all_zeros, 1, 2.8);
	}
	if (powered && first_pass) clear_colors_and_disable_light(ENGINE_DEF_LIGHT);
}


void uobj_draw_data::draw_dwexterm() const {

	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv3(get_ndiv(ndiv/3)), ndiv4(get_ndiv(ndiv/4));
	setup_draw_ship();
	colorRGBA const tex_color(0.3, 0.4, 0.3, 1.0);
	uniform_scale(1.5);
	glTranslatef(0.0, 0.0, -0.15);

	// backbone + body
	select_texture(BCUBE2_TEX);
	draw_cube(point(0.0, 0.01,  0.45), 0.12, 0.14, 1.5, 1, 1, 2.0, 1);
	draw_cube(point(0.0, 0.01, -1.15), 0.12, 0.14, 0.1, 1, 1, 2.0, 1);
	draw_cube(point(0.0, 0.00, -0.70), 0.36, 0.16, 0.8, 1, 1, 2.0, 1);

	// the bridge
	set_ship_texture(SHIP_HULL_TEX);
	glPushMatrix();
	glTranslatef(0.0, 0.12, 1.2);
	glScalef(1.0, 1.0, 1.8);
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	draw_cylinder(0.02, 0.055, 0.04, ndiv2, 1);
	glPopMatrix();
	end_ship_texture();
	LT_GRAY.do_glColor();

	for (unsigned i = 0; i < 6; ++i) { // "ribs"
		glPushMatrix();
		glTranslatef(0.0, -0.18, 0.04+0.1*i);
		glScalef(1.4, 0.85, 1.0);
		draw_torus(0.02, 0.25, ndiv4, ndiv);
		glPopMatrix();
	}
	glPushMatrix();
	glTranslatef(0.0, -0.38, 0.01);
	glScalef(1.2, 0.8, 1.0);
	draw_cylinder(0.56, 0.028, 0.028, ndiv2, 1);
	glPopMatrix();
	tex_color.do_glColor();
	set_ship_texture(BCUBE2_TEX);

	// forward deck
	draw_cube(point(0.0, 0.11, 0.30), 0.32, 0.06, 0.60, 1, 1, 1.0, 1);

	for (unsigned i = 0; i < 2; ++i) { // forward "wings"
		float const val(i ? 1.0 : -1.0);
		glPushMatrix();
		glTranslatef(-0.279*val, -0.005, -0.01);
		glRotatef(-43.0*val, 0.0, 0.0, 1.0);
		glScalef(0.15, 1.0, 1.0);
		glRotatef(21.0*val, 0.0, 0.0, 1.0);
		draw_cylinder(0.62, 0.2, 0.2, 6, 1);
		glPopMatrix();
	}

	// rear deck
	draw_cube(point(0.0, 0.11, -0.70), 0.40, 0.06, 0.80, 1, 1, 1.0, 1);
	glPushMatrix();
	glTranslatef(0.0, 0.1, -1.2);
	glScalef(1.0, 0.13, 1.0);
	draw_cylinder(0.1, 0.1, 0.22, 6, 1);
	glPopMatrix();

	for (unsigned i = 0; i < 2; ++i) { // rear "wings"
		float const val(i ? 1.0 : -1.0);
		glPushMatrix();
		glTranslatef(-0.371*val, -0.082, -1.04);
		glRotatef(-43.0*val, 0.0, 0.0, 1.0);
		glRotatef(-13.5, 1.0, 0.0, 0.0);
		glScalef(0.08, 1.0, 1.0);
		glRotatef(24.0*val, 0.0, 0.0, 1.0);
		draw_cylinder(0.78, 0.3, 0.1, 6, 1);
		glPopMatrix();
	}

	// front
	glPushMatrix();
	glTranslatef(0.0, 0.052, 0.75);
	glScalef(1.0, 0.7, 1.0);
	draw_cylinder(0.4, 0.0, 0.1, ndiv, 1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 0.04, 1.2);
	glScalef(1.0, 1.0, 1.8);
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	draw_cylinder(0.08, 0.14, 0.06, ndiv, 1);
	glPopMatrix();
	color_a.do_glColor();
	end_ship_texture();

	// forward weapons
	glPushMatrix();
	glTranslatef(0.0, 0.0, 1.2);
	draw_cylinder(0.08, 0.04, 0.04, ndiv2, 1);
	glPopMatrix();

	if (ndiv > 8) { // deck weapons
		DK_GRAY.do_glColor(); // front, 2x forward, 2x top, rear, bottom, head
		point const locs[8] = {point(0.12,  0.14,  0.50), point(0.10,  0.14,  0.28), point(0.10,  0.14,  0.12),
			                   point(0.15,  0.14, -0.55), point(0.15,  0.14, -0.85), point(0.55, -0.29, -0.97),
							   point(0.00, -0.41,  0.29), point(0.00,  0.14,  1.2)};
		for (unsigned i = 0; i < 8; ++i) {
			unsigned const nx(1 + (locs[i].x != 0.0));
			if (i == 6) LT_GRAY.do_glColor();

			for (unsigned j = 0; j < nx; ++j) {
				bool const half(i != 5 && i != 6);
				float const size((i==0) ? 0.042 : (half ? 0.030 : 0.024));
				glPushMatrix();
				point p(locs[i]);
				if (j) p.x = -p.x;
				translate_to(p);
				if (half) glRotatef(-90.0, 1.0, 0.0, 0.0);
				if (half) glScalef(1.0, 1.0, 0.6);
				draw_sphere_vbo(all_zeros, size, ndiv3, 0, half);
				glPopMatrix();
			}
		}
	}
	if (ndiv > 4) {
		if (ndiv > 6) {
			LT_GRAY.do_glColor();

			for (unsigned i = 0; i < 8; ++i) { // weapon barrels
				float const theta(TWO_PI*i/8.0), x(0.028*sinf(theta)), y(0.028*cosf(theta)); // are x and y backwards?
				glPushMatrix();
				glTranslatef(x, y, 1.28);
				draw_cylinder(0.06, 0.006, 0.006, ndiv4, (ndiv > 8));
				glPopMatrix();
			}
			for (unsigned i = 0; i < 2; ++i) { // connectors
				draw_cube(point((i ? 0.18 : -0.18), 0.10, -0.15), 0.02, 0.02, 0.30, 0);
			}
		}
		GRAY.do_glColor();
		set_ship_texture(SHIP_HULL_TEX);

		for (unsigned i = 0; i < 4; ++i) { // engines
			draw_ehousing_pairs(0.6-0.08*i, 0.04, 0.0, 0.034, 0.48+0.14*i, 0.0, 1, point(-0.24-0.07*i, 0.05-0.07*i, -1.12+0.02*i));
		}
		end_ship_texture();
	}
	if (first_pass) {
		assert(obj);
		unsigned const ndiv23(get_ndiv((2*ndiv)/3));
		float const instability(0.012 + ((t_exp == 0.0) ? 0.0 : 0.5*(1.0 - t_exp)));
		point lpos0(point(0.0, -0.205, -0.14)*1.5), lpos(lpos0);
		obj->xform_point_inv(lpos);
		glTranslatef(0.0, -0.205, 0.29);
		uniform_scale(0.5);
		draw_usw_star_int(ndiv23, lpos, lpos0, 2.0, 0.25, instability, powered); // will do a glPopMatrix()

		if (t_exp == 1.0 && animate2) { // exploding
			add_blastr(lpos, dir, 2.0*radius, 0.0, int(1.5*TICKS_PER_SECOND), -1, WHITE, WHITE, ETYPE_FUSION, obj);
		}
	}
	else {
		glPopMatrix(); // undo invert_z()
	}

	// draw engine glow
	draw_engine_pairs(BLUE, 0, 0.24, 0.36, 0.07, 1.96, point(0.1, -0.1, -0.03), 4);
}


void uobj_draw_data::draw_wraith_tail(float r, int ndiv2, float rscale) const {

	rscale = 0.75*fabs(rscale) + 0.25;
	point last_pos;

	for (unsigned i = 0; i < 26; ++i) { // default is a semicircle in the yz plane
		float const val(i/25.0), iscale(val*rscale + (1.0 - val)), iscal_inv(1.2/iscale);
		float const theta(iscale*TWO_PI*(i + 4.0)/36.0), rs(r*(1.0 - 0.02*i));
		point const pos(0.0, iscal_inv*(sinf(theta)-0.5), iscal_inv*(cosf(theta)-0.75));

		if (i == 25) {
			vector3d const dir(pos, last_pos);
			draw_fast_cylinder(last_pos, (pos + dir*4.0), rs, 0.0, ndiv2, 1);
		}
		else {
			if (ndiv <= 12) {
				draw_sphere_vbo(pos, rs, ndiv2, 1);
			}
			else {
				draw_sphere_vbo_back_to_front(pos, rs, ndiv2, 1);
			}
		}
		last_pos = pos;
	}
}


void uobj_draw_data::draw_wraith() const { // use time and vel_orient, fix bounding volume

	unsigned const ndiv_head(get_ndiv(3*ndiv/4)), ndiv_tail(get_ndiv(ndiv/3));
	setup_draw_ship();

	if (phase1) { // draw head - eyes?
		glPushMatrix();
		glScalef(0.8, 0.6, 1.0);
		draw_sphere_vbo(point(0.0, 0.22, 0.9), 0.6, ndiv_head, 0);
		glPopMatrix();
	}

	// draw body
	select_texture(NOISE_TEX);
	colorRGBA color_c(color_b);
	color_c.alpha = 0.8;
	color_c.do_glColor(); // bluegreen

	if (phase1) {
		glPushMatrix();
		glScalef(1.3, 0.5, 1.0);

		if (ndiv <= 10) {
			draw_sphere_vbo(point(0.0, 0.0, 0.75), 0.75, ndiv, 1);
		}
		else {
			draw_sphere_vbo_back_to_front(point(0.0, 0.0, 0.75), 0.75, ndiv, 1);
		}
		glPopMatrix();
	}
	//vector3d vel_orient(0.0, 0.0, 1.0);
	//if (vel != zero_vector) rotate_vector3d_by_vr(dir, vel.get_norm(), vel_orient);
	
	// draw tail
	float const tail_speed(0.005), rscale(cosf(tail_speed*TWO_PI*(on_time%unsigned(1.0/tail_speed)))); // -1.0 to 1.0

	if (ndiv > 6) {
		if (phase1) draw_wraith_tail(0.19, ndiv_tail, rscale);
		color_b.do_glColor();
	}
	if (phase2) draw_wraith_tail(0.27, ndiv_tail, rscale);
	end_texture();
	glPopMatrix();
	if (is_moving() && phase1) add_light_source(pos, 4.0*radius, color_b); // radius = f(health)?
}


void uobj_draw_data::draw_abomination() const {

	assert(obj);
	float const val(obj->get_state_val()); // 0.0 => fully closed, 1.0 => fully open
	unsigned const ndiv32(3*ndiv/2), eyelid_ndiv(min(32, 2*ndiv));
	setup_draw_ship();
	if (specular_en) set_specular(0.8, 80.0);
	glTranslatef(0.0, 0.0, 2.0);

	if (val > 0.0) {
		WHITE.do_glColor();
		draw_sphere_vbo(point(0.0, 0.0, 0.45), 0.8, ndiv, 0); // eyeball

		// rotate pupil to face the target when the eye is open and the target is in front
		// Note: not quite correct, since direction is measured from abomination center, not eyeball
		color_a.do_glColor();
		vector3d const tdir_inv_z(-tdir.x, tdir.y, -tdir.z); // invert_z() inverts x and z
		vector3d const pupil_dir((dot_product(tdir_inv_z, plus_z) > 0.2) ? tdir_inv_z : plus_z);
		draw_sphere_vbo((point(0.0, 0.0, 0.45) + 0.5*pupil_dir), 0.4, get_ndiv(3*ndiv/4), 0, 0); // pupil
	}
	color_b.do_glColor();

	if (0) {
		draw_subdiv_sphere_section(point(0.0, 0.0, 0.5), 1.0, eyelid_ndiv, 0, 0.0, 1.0, 0.4*val, 1.0); // eye hole
	}
	else {
		glPushMatrix();
		glRotatef(-90.0, 0.0, 1.0, 0.0);
		glRotatef(-90.0, 0.0, 0.0, 1.0);
		draw_subdiv_sphere_section(point(0.0, 0.5, 0.0), 1.0, eyelid_ndiv, 0, 0.18*val, (1.0-0.18*val), 0.0, 1.0); // eye slit
		glPopMatrix();
	}
	glPopMatrix(); // undo transformations
	assert(obj);
	cobj_vector_t const &cobjs(obj->get_cobjs());
	assert(!cobjs.empty() || on_time == 0);

	for (unsigned i = 2; i < cobjs.size(); ++i) { // skip the first one (the main sphere)
		cobjs[i]->draw(max(3, ((i == 2) ? 2*ndiv : (int(ndiv32) - int(i)))));
	}
	end_specular();
	
	if (powered && val > 0.1) {
		colorRGBA color;
		blend_color(color, WHITE, color_a, 0.5, 1);
		add_light_source((pos + dir*(3.5*radius)), 3.5*val*radius, color); // eye light
	}
}


void uobj_draw_data::draw_reaper() const {

	if (specular_en) set_specular(0.9, 90.0);
	setup_draw_ship();
	if (can_have_engine_lights()) setup_point_light(all_zeros, color_b, 5.0*radius, ENGINE_DEF_LIGHT);
	cobj_vector_t const &cobjs(obj->get_cobjs());

	if (cobjs.size() == 2) { // blocking shield is up
		color_a.do_glColor(); // never cloaked (semi-transparent? have to deal with draw order)
		glPushMatrix();
		invert_z(); // invert back to original orientation
		cobjs[1]->draw(ndiv);
		glPopMatrix();
	}
	color_b.do_glColor();
	//set_ship_texture(NOISE_TEX);
	draw_sphere_vbo_back_to_front(all_zeros, 1.0, 3*ndiv/2, 0);
	//end_ship_texture();
	if (can_have_engine_lights()) clear_colors_and_disable_light(ENGINE_DEF_LIGHT);
	end_specular();
	glPopMatrix(); // undo invert_z()
	glPopMatrix(); // undo rotations
	if (powered) setup_colors_draw_flare(pos, all_zeros, 3.0, 3.0, color_b);
	glPushMatrix();
	//if (powered && animate2 && final_pass && phase1) add_colored_lights(pos, radius, color_a, 0.25, 4, obj);
}


void uobj_draw_data::draw_death_orb() const {

	if (specular_en) set_specular(0.9, 90.0);
	setup_draw_ship();
	if (ndiv > 3) draw_sphere_vbo(all_zeros, 0.5, ndiv, 0);
	end_specular();
	glPopMatrix(); // undo transformations

	if (powered) {
		setup_colors_draw_flare(pos, all_zeros, 2.1, 2.1, color_a, FLARE2_TEX);
		if (ndiv > 3) add_light_source(pos, 6.0*radius, color_a);
	}
	if (is_moving() && animate2 && first_pass && ndiv > 4) {
		bool const two_parts((time%5)==0);

		for (unsigned i = 0; i < unsigned(1+two_parts); ++i) {
			gen_particle(PTYPE_GLOW, color_a, ALPHA0, TICKS_PER_SECOND, pos, zero_vector,
				0.15*radius, 0.0, obj->get_align(), (i > 0));
		}
	}
}


void uobj_draw_data::draw_supply() const {

	unsigned const ndiv2(get_ndiv(ndiv/2));
	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	draw_cylinder(point(0.0, 0.0, 1.0), 0.4, 0.5, 0.5, ndiv, 1); // front ring
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.0, -1.0);
	glScalef(1.0, 1.8, 1.0);
	draw_cylin_fast(0.1, 0.1, 1.8, ndiv2, textured); // backbone
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.0, 0.0, 0.8);
	draw_cylin_fast(0.45, 0.45, 0.8, ndiv, textured); // front cylinder
	glTranslatef(0.0, 0.0, -0.6);
	draw_cylin_fast(0.00, 0.45, 0.6, ndiv, textured);
	glTranslatef(0.0, 0.0, -1.2);
	glScalef(1.0, 1.5, 1.0);
	draw_cylin_fast(0.3, 0.0, 0.6, ndiv, textured);
	glPopMatrix();
	glPushMatrix();
	glScalef(1.0, 1.0, 0.5);
	draw_sphere_vbo(point(0.0, 0.0, 3.2), 0.45, ndiv, textured, 1); // front sphere
	glPopMatrix();
	draw_cube(point(0.0, 0.0, -1.3), 0.6, 0.9, 0.6, textured, 1); // rear

	if (ndiv > 4) { // draw engines
		draw_ehousing_pairs(0.55, 0.2, 0.12, 0.1, 0.5, 0.0, 1, point(-0.25, -0.4, -1.7), point(0.0, 0.4, 0.0), 3);
	}
	if (textured) end_ship_texture();
	colorRGBA light_color(BLACK);

	if (powered && ((2*time/TICKS_PER_SECOND) & 1)) { // draw blinky light
		for (unsigned i = 0; i < 3; ++i) {
			light_color[i] = ((color_a[i] > 0.8) ? 1.0 : 0.0);
		}
		if (light_color == BLACK) light_color = WHITE;
	}
	set_emissive_color(light_color);
	draw_sphere_vbo(point(0.0, 0.0, 1.825), 0.07, get_ndiv(ndiv/4), 0);
	clear_emissive_color();
	glPopMatrix(); // undo invert_z()

	for (unsigned i = 0; i < 3; ++i) {
		draw_engine_pairs(LT_BLUE, (i<<1), 0.5, 0.25, (-0.4 + 0.4*i), 1.75);
	}
}


void uobj_draw_data::draw_anti_miss() const {

	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv3(get_ndiv(ndiv/3));
	setup_draw_ship();
	glTranslatef(0.0, 0.36, 0.0);
	draw_sphere_vbo(all_zeros, 0.75, 3*ndiv/2, 0); // main body
	color_b.do_glColor();
	
	if (ndiv > 6) { // weapon
		glPushMatrix();
		invert_z();
		rotate_into_plus_z(tdir);
		draw_cylin_fast(0.08, 0.06, 1.2, ndiv3, 0);
		glPopMatrix();
	}
	glPushMatrix();
	glTranslatef(0.0, -0.5, 0.0);
	glRotatef(90.0, 1.0, 0.0, 0.0);
	draw_torus(0.2, 0.9, ndiv2, 4*ndiv/3); // ring
	glPopMatrix();

	if (ndiv > 4) {
		for (unsigned i = 0; i < 3; ++i) { // tripod
			float const theta(i*TWO_PI/3.0);
			point const pt(0.84*cosf(theta), -0.8, 0.84*sinf(theta));
			color_a.do_glColor();
			draw_sphere_vbo(pt, 0.25, ndiv2, 0);
			
			if (ndiv > 8) {
				color_b.do_glColor();
				draw_fast_cylinder(point(0.0, -0.7, 0.0), pt, 0.07, 0.07, ndiv3, 0, 0);
			}
		}
	}
	glPopMatrix();
}


void uobj_draw_data::draw_juggernaut() const {

	unsigned const ndiv2(get_ndiv(ndiv/2));
	setup_draw_ship();
	if (powered && first_pass) setup_point_light(point(0.0, 0.05, -0.95), color_b, 1.0*radius, ENGINE_DEF_LIGHT);

	glPushMatrix();
	glScalef(0.85, 1.2, 1.0);
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 0); // main body
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0, 0.5, 0.3);
	glScalef(0.5, 1.2, 1.2);
	draw_sphere_vbo(all_zeros, 0.6, ndiv, 0); // head
	glPopMatrix();

	if (ndiv > 8) {
		RED.do_glColor();
		unsigned const ndiv4(get_ndiv(ndiv/4));
		point eyes[3] = {point(0.0, 0.76, 0.96), point(-0.08, 0.6, 0.98), point(0.08, 0.6, 0.98)};

		for (unsigned i = 0; i < 3; ++i) { // eyes
			draw_sphere_vbo(eyes[i], 0.04, ndiv4, 0); // head
		}
	}
	glPushMatrix();
	glTranslatef(0.0, 0.05, -0.45);
	glScalef(0.55, 1.5, 1.0);
	if (powered) set_emissive_color(color_b); else color_b.do_glColor();
	draw_sphere_vbo(all_zeros, 0.65, ndiv, 0); // back
	if (powered) clear_emissive_color();
	glPopMatrix();
	color_a.do_glColor();

	for (unsigned i = 0; i < 2; ++i) {
		float const is(i ? 1.0 : -1.0);
		draw_sphere_vbo(point(0.64*is, -0.8, -0.25), 0.55, ndiv, 0); // bottom
		
		glPushMatrix();
		glTranslatef(0.9*is, -0.2, 0.5);
		glScalef(0.8, 0.5, 2.3);
		draw_sphere_vbo(all_zeros, 0.35, ndiv2, 0); // middle
		glPopMatrix();
		
		RED.do_glColor();
		glPushMatrix();
		glTranslatef(0.9*is, -0.2, 0.65);
		glScalef(0.7, 0.5, 2.8);
		draw_sphere_vbo(all_zeros, 0.28, ndiv2, 0); // red part
		glPopMatrix();

		color_a.do_glColor();
		glPushMatrix();
		glTranslatef(0.66*is, 0.7, 0.1);
		glScalef(0.9, 0.9, 1.6);
		draw_sphere_vbo(all_zeros, 0.5, ndiv, 0); // top
		glPopMatrix();
	}
	glPopMatrix();
	if (powered && first_pass) clear_colors_and_disable_light(ENGINE_DEF_LIGHT);
}


void uobj_draw_data::draw_saucer(bool rotated, bool mothership) const {

	unsigned const ndiv32(3*ndiv/2), ndiv2(get_ndiv(ndiv/2)), ndiv4(get_ndiv(ndiv/4));
	bool const pt_light(can_have_engine_lights() && !(eflags & 1));
	setup_draw_ship();
	uniform_scale(1.5);
	if (rotated) glRotatef(-90.0, 1.0, 0.0, 0.0); // rotate so that "up" is in +y
	if (mothership) glScalef(1.0, 1.0, 0.75);
	color_b.do_glColor(); // WHITE?
	if (specular_en) set_specular(0.9, 90.0);
	if (pt_light) setup_point_light(point(0.0, 0.0, -0.5), color_a, 3.0*radius, ENGINE_DEF_LIGHT);
	set_ship_texture(SPACESHIP2_TEX);
	glPushMatrix();
	glScalef(1.0, 1.0, 0.1);
	draw_sphere_vbo(all_zeros, 1.0, ndiv32, 1); // center
	glPopMatrix();
	draw_cylin_fast(0.95, 0.15, 0.4, ndiv32, 1); // top
	glPushMatrix();
	glTranslatef(0.0, 0.0, -0.3);
	draw_cylin_fast(0.2, 0.8, 0.3, ndiv32, 1); // bottom
	glPopMatrix();
	end_ship_texture();
	if (pt_light) clear_colors_and_disable_light(ENGINE_DEF_LIGHT);

	if (ndiv > 4) {
		color_a.do_glColor();
		glPushMatrix();
		glTranslatef(0.0, 0.0, 0.4);

		if (mothership) {
			draw_cylin_fast(0.15, 0.0, 0.6, ndiv2, 0); // pointy top cylinder
		}
		else {
			glScalef(1.0, 1.0, 0.5);
			draw_sphere_vbo(all_zeros, 0.15, ndiv2, 0); // topmost sphere
		}
		glPopMatrix();
	}
	if (ndiv > 6) {
		unsigned const nwpts(mothership ? 10 : 6);

		for (unsigned i = 0; i < nwpts; ++i) {
			((i&1) ? colorRGBA(0.5, 0.3, 0.3, 1.0) : (mothership ? colorRGBA(0.8, 0.5, 0.2, 1.0) : colorRGBA(0.3, 0.3, 0.5, 1.0))).do_glColor();
			float const theta(TWO_PI*i/float(nwpts));

			draw_sphere_vbo(point(0.6*cosf(theta), 0.6*sinf(theta), 0.16), 0.065, ndiv2, 0);
		}
	}
	if (ndiv > 8) {
		unsigned const nlights(mothership ? 16 : 12);

		for (unsigned is_lit = 0; is_lit < 2; ++is_lit) {
			if (is_lit) {
				set_emissive_color(RED); // always red
			}
			else {
				color_b.do_glColor();
			}
			for (unsigned i = 0; i < nlights; ++i) {
				if ((powered && ((i+(on_time>>2))&3) == 0) != is_lit) continue; // incorrect state
				float const theta(TWO_PI*i/float(nlights));
				draw_sphere_vbo(point(cosf(theta), sinf(theta), 0.0), 0.045, ndiv4, 0);
			}
			if (is_lit) clear_emissive_color();
		}
	}
	end_specular();
	glPopMatrix(); // undo invert_z()

	if (ndiv > 4 && is_moving() && !(eflags & 1) && vel.mag_sq() > 1.0E-7) { // draw engine glow
		enable_ship_flares(color_a);
		draw_engine(color_a, (rotated ? point(0.0, -0.5, 0.0) : point(0.0, 0.0, 0.5)), 1.0);
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_headhunter() const {

	setup_draw_ship();

	if (ndiv > 3) { // draw fins
		vert_norm verts[12];

		for (unsigned dim = 0; dim < 2; ++dim) { // {x,y}
			for (unsigned dir = 0; dir < 2; ++dir) { // {-,+}
				float const val(dir ? 1.0 : -1.0);
				vector3d const n(val*dim, val*(1-dim), 0);
				verts[3*(2*dim+dir)+0] = vert_norm(point(0.05*val*(1-dim), 0.05*val*dim, -0.6), n);
				verts[3*(2*dim+dir)+1] = vert_norm(point(0.05*val*(1-dim), 0.05*val*dim, -1.6), n);
				verts[3*(2*dim+dir)+2] = vert_norm(point(0.45*val*(1-dim), 0.45*val*dim, -1.7), n);
			}
		}
		draw_verts(verts, 12, GL_TRIANGLES);
	}

	// draw body
	if (specular_en) set_specular(0.9, 90.0);
	color_b.do_glColor();
	glScalef(0.2, 0.2, 1.7); // Note: push/pop not needed since this is the last draw
	draw_sphere_vbo(point(0.0, 0.0, 0.03), 1.0, ndiv, 0);
	if (specular_en) set_specular(0.0, 0.0);
	glPopMatrix(); // undo invert_z()

	if (is_moving() && !(eflags & 1)) { // draw engine glow
		colorRGBA const color(0.75, 0.75, 1.0, 1.0);
		enable_ship_flares(color);
		draw_engine(color, point(0.0, 0.0, 1.8), 0.8);
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_seige() const {

	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv4(get_ndiv(ndiv/2));
	setup_draw_ship();
	select_texture(NOISE_TEX);

	// draw top
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.0, 0.4);
	glScalef(0.6, 0.2, 1.0);
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 1);
	glTranslatef(0.0, 0.0, -0.3);
	glScalef(1.0, 1.0, 1.5);
	colorRGBA c(color_b);

	for (unsigned i = 0; i < 4; ++i) {
		c *= 0.9;
		c.do_glColor();
		glScalef(0.7, 1.35, 0.8);
		glTranslatef(0.0, 0.05*(7.0 - i), 0.0);
		draw_sphere_vbo(all_zeros, 1.0, ndiv, 1);
	}
	glPopMatrix();

	// draw rear
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.0, -0.4);
	glScalef(0.9, 0.25, 0.96);
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 1);
	glPopMatrix();

	// draw engines
	point epts[3];
	glTranslatef(0.0, 0.0, -0.4);

	for (unsigned i = 0; i < 3; ++i) { // create engines
		float const theta((i - 1.0)*PI/5.0);
		epts[i] = point(sinf(theta), 0.0, -cosf(theta));
		if (ndiv < 8) continue;
		glPushMatrix();
		translate_to(epts[i]);
		glRotatef(TO_DEG*theta, 0.0, -1.0, 0.0);
		glScalef(1.0, 0.4, 1.0);
		draw_cylin_fast(0.2, 0.2, 0.3, ndiv2, 1); // engine housing
		glPopMatrix();
	}

	// draw bottom
	color_b.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.1, 0.4);
	glScalef(0.3, 1.0, 1.0);
	glRotatef(70.0, 1.0, 0.0, 0.0);
	draw_cylin_fast(0.7, 0.0, 1.3, ndiv, 1);
	end_texture();
	color_a.do_glColor();
	glScalef(1.0/0.3, 1.0, 1.0);
	draw_sphere_vbo(point(0.0, 0.0, 1.3), 0.12, ndiv2, 0); // weapon sphere
	glPopMatrix();

	// draw weapons
	if (ndiv > 3) {
		unsigned wpt(0);
		float ftime(-1.0);

		if (ndiv >= 10 && obj != NULL) { // animate + flash
			vector<ship_weapon> const &weapons(*(obj->get_weapons()));

			for (unsigned i = 0; i < weapons.size(); ++i) {
				if (weapons[i].wclass == UWEAP_SEIGEC) {
					wpt   = weapons[i].cur_wpt;
					ftime = float(frame_counter - weapons[i].last_fframe)/TICKS_PER_SECOND;
					break;
				}
			}
		}
		for (unsigned i = 0; i < 2; ++i) {
			float comp((wpt == i && ftime >= 0.0 && ftime < 1.0) ? (0.5 + 0.5*ftime) : 1.0);
			glPushMatrix();
			glTranslatef(0.13*(i ? -1.0 : 1.0), 0.6, 0.5);
			draw_cylinder(0.32, 0.032, 0.032, ndiv4, 1);

			for (unsigned j = 0; j < 4; ++j) {
				glTranslatef(0.0, 0.0, 0.24*comp);
				float const r(0.028 - 0.004*j);
				draw_cylin_fast(r, r, 0.32, ndiv4, 0);
			}
			glPopMatrix();
		}
	}
	glPopMatrix(); // undo transformations

	if (is_moving()) { // draw engine glow
		colorRGBA const color(blend_color(color_a, WHITE, 0.5, 1));
		enable_ship_flares(color);

		for (unsigned i = 0; i < 3; ++i) {
			if (eflags & (1 << i)) continue;
			vector3d edir(cross_product(epts[i], vector3d(0.0, 1.0, 0.0)).get_norm());
			point ept(epts[i]*1.05);
			edir.z *= -1.0;
			ept.z   = -ept.z + 0.4;
			draw_engine(color, ept, 0.3, 2.5, edir);
		}
		disable_ship_flares();
	}
}


void uobj_draw_data::draw_colony(bool armed, bool hw, bool starport) const {

	unsigned const ndiv2(get_ndiv(ndiv/2)), ndiv4(get_ndiv(ndiv/4));
	setup_draw_ship();
	bool const textured(1);
	if (textured) set_ship_texture(SHIP_HULL_TEX);

	glPushMatrix();
	glTranslatef(0.0, 0.0, -1.0);
	draw_cylinder(1.5, 1.0, 1.0, 3*ndiv/2, 1, 1, 0);
	glPopMatrix();

	glPushMatrix();
	glScalef(1.0, 1.0, 0.5);
	draw_sphere_vbo(point(0.0, 0.0, 1.0), 1.0, 3*ndiv/2, textured);
	glPopMatrix();
	color_b.do_glColor();

	if (armed && ndiv > 4) { // draw weapon
		glPushMatrix();
		glTranslatef(0.0, 0.0, 1.0);
		draw_sphere_vbo(all_zeros, 0.15, ndiv2, 0);

		if (ndiv > 5) {
			invert_z();
			vector3d turret_dir(tdir);
			turret_dir.z = min(0.0f, turret_dir.z);
			rotate_into_plus_z(turret_dir);
			draw_cylin_fast((hw ? 0.04 : 0.02), (hw ? 0.032 : 0.016), (hw ? 0.9 : 0.6), ndiv4, 0);
		}
		glPopMatrix();
	}
	if (hw && ndiv > 3) {
		for (unsigned i = 0; i < 3; ++i) {
			float const theta(TWO_PI*i/3.0), x(1.05*cosf(theta)), y(1.05*sinf(theta));
			glPushMatrix();
			glTranslatef(x, y, -0.6);
			draw_cylinder(1.1, 0.2, 0.2, ndiv2, 1, 1, 0);
			draw_sphere_vbo(point(0.0, 0.0, 0.0), 0.2, ndiv2, textured);
			draw_sphere_vbo(point(0.0, 0.0, 1.1), 0.2, ndiv2, textured);
			glPopMatrix();
		}
	}
	if (starport && ndiv > 3) {
		glTranslatef(0.0, 0.0, 0.25); // Note: push/pop not needed since this is the last draw
		draw_torus(0.2, 1.05, ndiv2, 3*ndiv/2);
	}
	if (textured) end_ship_texture();
	glPopMatrix();
}


void uobj_draw_data::draw_default_ship() const {

	setup_draw_ship();
	draw_sphere_vbo(all_zeros, 1.0, 2*ndiv, 0);
	glPopMatrix(); // undo transformations
}


// ******************* STELLAR OBJECTS *******************


void uobj_draw_data::draw_asteroid(int tex_id) const {

	color_a.do_glColor();
	select_texture(tex_id);
	draw_sphere_vbo(all_zeros, 1.0, 3*ndiv/2, 1);
	end_texture();
}


void uobj_draw_data::draw_black_hole() const { // should be non-rotated

	BLACK.do_glColor();
	draw_sphere_vbo(all_zeros, 0.3, ndiv, 0);
	vector3d const player_dir(player_ship().get_pos() - pos);
	float const dist_to_player(player_dir.mag());

	if (dist_to_player > 2.0*radius) {
		quad_batch_draw qbd;
		qbd.add_xlated_billboard(pos, player_dir*(1.2/dist_to_player), get_camera_pos(), up_vector, BLACK, 3.0, 3.0);
		qbd.draw_as_flares_and_clear(BLUR_TEX);
	}
}




