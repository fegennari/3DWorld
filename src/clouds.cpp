// 3D World - Cloud/Nebula Generation and Drawing Code
// by Frank Gennari
// 3/10/12
#include "3DWorld.h"
#include "physics_objects.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "draw_utils.h"
#include "clouds.h" // for uobj_base for unebula


bool     const USE_CLOUD_FBO    = 1;
unsigned const CLOUD_GEN_TEX_SZ = 1024;


vector2d cloud_wind_pos(0.0, 0.0);
cloud_manager_t cloud_manager;

extern bool have_sun, no_sun_lpos_update;
extern int window_width, window_height, cloud_model, draw_model, display_mode, xoff, yoff, animate2;
extern float CLOUD_CEILING, atmosphere, sun_rot, fticks, water_plane_z, zmin, zmax;
extern vector3d wind;


void draw_part_cloud(vector<particle_cloud> const &pc, colorRGBA const color, bool zoomed);


struct cloud_t {
	unsigned b, e;
	point p;
	float r;
	cloud_t(unsigned b_=0, unsigned e_=0, point const &p_=all_zeros, float r_=0.0)
		: b(b_), e(e_), p(p_), r(r_) {}
};


void cloud_manager_t::create_clouds() { // 3D cloud puffs

	if (!empty()) return; // keep the old clouds
	clear();
	free_textures();
	srand(123);
	float const xsz(X_SCENE_SIZE), ysz(Y_SCENE_SIZE);
	unsigned const NCLOUDS = 10;
	unsigned const NPARTS  = 1000;

	for (unsigned c = 0; c < NCLOUDS; ++c) {
		point const center(4.0*xsz*signed_rand_float(), 4.0*ysz*signed_rand_float(),
						   (ztop + CLOUD_CEILING + Z_SCENE_SIZE*rand_uniform(0.25, 0.75)));
		point const bounds(xsz*rand_uniform(1.0, 2.0), ysz*rand_uniform(1.0, 2.0),
						   Z_SCENE_SIZE*rand_uniform(0.4, 0.8));
		unsigned const nparts(rand()%(NPARTS/2) + NPARTS/2);
		size_t const ix(size());
		resize(ix + nparts);

		for (unsigned p = 0; p < nparts; ++p) {
			point pos(signed_rand_vector_spherical(1.0));

			for (unsigned i = 0; i < 3; ++i) {
				pos[i] *= bounds[i];
			}
			if (pos.z < 0.0) pos.z *= 0.5; // compressed on the bottom
			pos += center;
			float const radius(0.045*(xsz + ysz)*rand_uniform(0.5, 1.0));
			float const density(rand_uniform(0.05, 0.12));
			(*this)[ix + p].gen(pos, WHITE, zero_vector, radius, density, 0.0, 0.0, -((int)c+2), 0, 0, 1, 1); // no lighting
		}
	}
}


void cloud_manager_t::update_lighting() {

	RESET_TIME;
	point const sun_pos(get_sun_pos());
	bool const calc_sun_light(have_sun && light_factor > 0.4);
	unsigned const num_clouds((unsigned)size());
	int last_src(0);
	vector<cloud_t> clouds;

	if (calc_sun_light) {
		for (unsigned i = 0; i < num_clouds; ++i) {
			particle_cloud &c((*this)[i]);

			if (i == 0 || c.source != last_src) {
				last_src = c.source;
				if (i > 0) clouds.back().e = i; // end the last cloud
				clouds.push_back(cloud_t(i));   // begin a new cloud
			}
		}
		clouds.back().e = num_clouds; // end the last cloud

		for (unsigned i = 0; i < clouds.size(); ++i) {
			cloud_t &c(clouds[i]);
			for (unsigned j = c.b; j < c.e; ++j) c.p += (*this)[j].pos;
			c.p /= (c.e - c.b);
			for (unsigned j = c.b; j < c.e; ++j) c.r = max(c.r, (p2p_dist(c.p, (*this)[j].pos) + (*this)[j].radius));
		}
	}
	for (unsigned i = 0; i < num_clouds; ++i) {
		particle_cloud &pc((*this)[i]);
		float light(0.25); // night time sky

		if (calc_sun_light) {
			vector3d const v1(sun_pos - pc.pos);
			float const dist_sq(v1.mag_sq());
			vector3d const v1n(v1/dist_sq);
			light = 1.0; // start off fully lit

			for (unsigned p = 0; p < clouds.size(); ++p) {
				cloud_t &c(clouds[p]);
				float t; // unused

				if (sphere_test_comp(sun_pos, c.p, v1, c.r*c.r, t)) {
					for (unsigned j = c.b; j < c.e; ++j) {
						particle_cloud const &c2((*this)[j]);
						vector3d const v2(sun_pos, c2.pos);
						if (v2.mag_sq() > dist_sq) continue; // further from the sun
						float const dotp(dot_product(v1, v2));
						float const dsq((dotp > dist_sq) ? p2p_dist_sq(v1, v2) : (v2 - v1n*dotp).mag_sq());
						if (dsq > c2.radius*c2.radius) continue; // no intersection
						float const alpha(2.0*c2.base_color.alpha*c2.density*((c2.radius - sqrt(dsq))/c2.radius));
						light *= 1.0 - CLIP_TO_01(alpha);
					}
				}
			}
			if (light_factor < 0.6) {
				float const blend(5.0*(light_factor - 0.4));
				light = light*blend + 0.25*(1.0 - blend);
			}
		}
		pc.darkness   = 1.0 - 2.0*light;
		pc.base_color = WHITE;
		apply_red_sky(pc.base_color);
	}
	PRINT_TIME("Cloud Lighting");
}


cube_t cloud_manager_t::get_bcube() const {

	cube_t bcube;

	for (unsigned i = 0; i < size(); ++i) {
		point const &pos((*this)[i].pos);
		float const radius((*this)[i].radius);

		if (i == 0) {
			bcube = cube_t(pos, pos);
			bcube.expand_by(radius);
		}
		else {
			bcube.union_with_sphere(pos, radius);
		}
	}
	return bcube;
}


float cloud_manager_t::get_max_xy_extent() const {

	cube_t const bcube(get_bcube());
	return(max(max(-bcube.d[0][0], bcube.d[0][1]), max(-bcube.d[1][0], bcube.d[1][1])));
}


bool cloud_manager_t::create_texture(bool force_recreate) {

	RESET_TIME;
	unsigned const xsize(USE_CLOUD_FBO ? CLOUD_GEN_TEX_SZ : min(CLOUD_GEN_TEX_SZ, (unsigned)window_width));
	unsigned const ysize(USE_CLOUD_FBO ? CLOUD_GEN_TEX_SZ : min(CLOUD_GEN_TEX_SZ, (unsigned)window_height));

	if (txsize != xsize || tysize != ysize) {
		free_textures();
		txsize = xsize;
		tysize = ysize;
	}
	if (cloud_tid && !force_recreate) return 0; // nothing to do
	
	if (!cloud_tid) {
		setup_texture(cloud_tid, GL_MODULATE, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, xsize, ysize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	}
	assert(glIsTexture(cloud_tid));
	check_gl_error(800);
	if (USE_CLOUD_FBO) enable_fbo(fbo_id, cloud_tid, 0);
	check_gl_error(801);

	glViewport(0, 0, xsize, ysize);
	glClearColor(1.0, 1.0, 1.0, 1.0); // white
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	// setup projection matrix
	cube_t const bcube(get_bcube());
	float const cloud_bot(bcube.d[2][0]), cloud_top(bcube.d[2][1]), cloud_xy(get_max_xy_extent());
	float const scene_xy(max(X_SCENE_SIZE, Y_SCENE_SIZE)), angle(atan2(cloud_xy, cloud_bot)), z1(min(zbottom, czmin));
	frustum_z = z1 - scene_xy*(cloud_bot - z1)/(cloud_xy - scene_xy);
	//pos_dir_up const pdu(get_pt_cube_frustum_pdu(get_camera_pos(), bcube, 1));
	//pos_dir_up const pdu(all_zeros, plus_z, plus_x, tanf(angle)*SQRT2, sinf(angle), NEAR_CLIP, FAR_CLIP, 1.0);
	//gluPerspective(2.0*angle/TO_RADIANS, 1.0, cloud_bot-frustum_z, cloud_top+(cloud_top - cloud_bot)-frustum_z);
	gluPerspective(2.0*angle/TO_RADIANS, 1.0, NEAR_CLIP, FAR_CLIP);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	vector3d const up_dir(plus_y);
	point const origin(0.0, 0.0, frustum_z), center(0.0, 0.0, cloud_bot);
	gluLookAt(origin.x, origin.y, origin.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);

	set_red_only(1);
	point const orig_cpos(camera_pos);
	bool const was_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable view frustum culling
	camera_pos = origin;
	draw_part_cloud(*this, WHITE, 1); // draw clouds
	camera_pos = orig_cpos;
	camera_pdu.valid = was_valid;
	set_red_only(0);

	if (!USE_CLOUD_FBO) { // render clouds to texture
		glBindTexture(GL_TEXTURE_2D, cloud_tid);
		glReadBuffer(GL_BACK);
		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, xsize, ysize); // copy the frame buffer to the bound texture
	}

	// reset state
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	if (USE_CLOUD_FBO) disable_fbo();
	glViewport(0, 0, window_width, window_height);
	if (!USE_CLOUD_FBO) glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	check_gl_error(802);
	PRINT_TIME("Cloud Texture Gen");
	return 1;
}


void free_cloud_textures() {

	cloud_manager.free_textures();
}


void cloud_manager_t::free_textures() {

	free_texture(cloud_tid);
	free_fbo(fbo_id);
}


//http://www.gamedev.net/reference/articles/article2273.asp
void cloud_manager_t::draw() {

	if (atmosphere < 0.01) return; // no atmosphere
	create_clouds();
	if (empty()) return;
	//glFinish(); // testing
	RESET_TIME;
	glDisable(GL_DEPTH_TEST);

	// WRITE: wind moves clouds

	// light source code
	static bool had_sun(0);
	static float last_sun_rot(0.0);
	bool const need_update(!no_sun_lpos_update && (sun_rot != last_sun_rot || have_sun != had_sun));

	if (need_update) {
		last_sun_rot = sun_rot;
		had_sun      = have_sun;
		update_lighting();
	}
	if (cloud_model == 0) { // faster billboard texture mode
		create_texture(need_update);
		enable_flares(get_cloud_color(), 1); // texture will be overriden
		assert(cloud_tid);
		bind_2d_texture(cloud_tid);

		shader_t s;
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("cloud_billboard");
		s.begin_shader();
		s.add_uniform_int("tex0", 0);

		glBegin(GL_QUADS);
		point const camera(get_camera_pos());
		cube_t const bcube(get_bcube());
		float const cloud_bot(bcube.d[2][0]), cloud_top(bcube.d[2][1]), cloud_xy(get_max_xy_extent());
		float const xy_exp((cloud_top - frustum_z)/(cloud_bot - frustum_z));
		
		for (unsigned d = 0; d < 2; ++d) { // render the bottom face of bcube
			for (unsigned e = 0; e < 2; ++e) {
				glTexCoord2f(float(d^e^1), float(d));
				point(xy_exp*((d^e) ? cloud_xy : -cloud_xy)+camera.x, xy_exp*(d ? cloud_xy : -cloud_xy)+camera.y, cloud_top).do_glVertex();
			}
		}
		glEnd();
		s.end_shader();
		disable_flares();
	}
	else {
		draw_part_cloud(*this, get_cloud_color(), 1);
	}
	glEnable(GL_DEPTH_TEST);
	//glFinish(); // testing
	//PRINT_TIME("Clouds");
}


void draw_puffy_clouds(int order) {

	if (cloud_manager.is_inited() && (get_camera_pos().z > cloud_manager.get_z_plane()) != order) return;

	if (atmosphere < 0.01) {
		cloud_manager.clear();
	}
	else if (display_mode & 0x40) { // key 7
		cloud_manager.draw();
	}
}


float get_cloud_zmax() {return get_camera_pos().z + max(zmax, CLOUD_CEILING);}


void set_cloud_uniforms(shader_t &s, unsigned tu_id) {

	select_multitex(NOISE_GEN_TEX, tu_id, 0);
	s.add_uniform_int("cloud_noise_tex", tu_id);
	set_multitex(0);
	s.add_uniform_vector2d("dxy", cloud_wind_pos);
}


void draw_cloud_plane(bool reflection_pass) {

	unsigned const NUM_DIV = 32;
	float const cloud_rel_vel = 1.0; // relative cloud velocity compared to camera velocity
	float const size(FAR_CLIP), rval(0.94*size), rval_inv(1.0/rval); // extends to at least the far clipping plane
	float const z1(zmin), z2(get_cloud_zmax()), ndiv_inv(1.0/NUM_DIV);
	point const camera(get_camera_pos()), world_pos(camera + vector3d((xoff2-xoff)*DX_VAL, (yoff2-yoff)*DY_VAL, 0.0));
	vector3d const offset(-camera + cloud_rel_vel*world_pos);
	static indexed_mesh_draw<vert_wrap_t> imd;
	shader_t s;
	glDepthMask(GL_FALSE);

	if (animate2) {
		cloud_wind_pos.x -= fticks*wind.x;
		cloud_wind_pos.y -= fticks*wind.y;
	}

	// draw a plane at zmin to properly blend the fog
	if (!reflection_pass) {
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		s.set_vert_shader("water_fog.part+fog_only");
		s.set_frag_shader("linear_fog.part+fog_only");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_float("water_plane_z", get_tiled_terrain_water_level());
		BLACK.do_glColor();
		imd.render_z_plane(-size, -size, size, size, zmin, NUM_DIV, NUM_DIV);
		s.end_shader();
	}

	// draw clouds
	s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
	s.set_prefix("#define NUM_OCTAVES 8",     1); // FS
	s.set_bool_prefix("underwater_atten", (glIsEnabled(GL_FOG) != 0), 1); // FS
	s.set_vert_shader("water_fog.part+clouds");
	s.set_frag_shader("linear_fog.part+perlin_clouds.part*+clouds");
	s.begin_shader();
	s.setup_fog_scale();
	s.add_uniform_float("water_plane_z", get_tiled_terrain_water_level());
	s.add_uniform_float("cloud_scale", 0.5);
	set_cloud_uniforms(s, 0);
	s.add_uniform_vector3d("cloud_offset", offset);
	enable_blend();
	get_cloud_color().do_glColor();
	imd.init(NUM_DIV, NUM_DIV);
	float xvals[NUM_DIV+1], yvals[NUM_DIV+1];

	for (unsigned t = 0; t <= NUM_DIV; ++t) {
		float const angle(TWO_PI*(t*ndiv_inv));
		xvals[t] = cosf(angle);
		yvals[t] = sinf(angle);
	}
	for (unsigned r = 0; r <= NUM_DIV; ++r) { // spherical section
		float const radius(size*r*ndiv_inv), zval(z1 + (z2 - z1)*cos(PI_TWO*min(1.0f, radius*rval_inv)));

		for (unsigned t = 0; t <= NUM_DIV; ++t) {
			imd.set_vert(t, r, point((radius*xvals[t] + camera.x), (radius*yvals[t] + camera.y), zval));
		}
	}
	imd.render();
	s.end_shader();
	disable_blend();
	glDepthMask(GL_TRUE);
}


// *** nebula code ***


unsigned const NUM_NEBULA_QUADS = 16;

void move_in_front_of_far_clip(point_d &pos, point const &camera, float &size, float dist);


void unebula::gen(float range, ellipsoid_t const &bounds) {

	// Note: bounds is not currently used, but it can be used to scale the nebula to the galaxy's ellipsoid (but requires some transforms in the shader)
	rand_gen_t rgen;
	rgen.set_state(rand2(), rand2());
	radius = rgen.rand_uniform(0.1, 0.15)*range;

	for (unsigned d = 0; d < 3; ++d) {
		color[d] = colorRGBA(rgen.rand_uniform(0.3, 1.0), rgen.rand_uniform(0.1, 0.5), rgen.rand_uniform(0.2, 0.9), 1.0);
	}
	points.resize(4*NUM_NEBULA_QUADS);

	for (unsigned i = 0; i < points.size(); i += 4) {
		vector3d const normal(rgen.signed_rand_vector_norm());
		vector3d vab[2];
		get_ortho_vectors(normal, vab);

		for (unsigned j = 0; j < 4; ++j) { // Note: quads will extend beyond radius, but will be rendered as alpha=0 outside radius
			points[i+j].v = radius*(((j>>1) ? 1.0 : -1.0)*vab[0] + (((j&1)^(j>>1)) ? 1.0 : -1.0)*vab[1]);
			points[i+j].set_norm(normal);
		}
	}
}


void unebula::begin_render(shader_t &s) {

	if (!s.is_setup()) {
		s.set_prefix("#define NUM_OCTAVES 5", 1); // FS
		s.set_bool_prefix("line_mode", (draw_model == 1), 1); // FS
		s.set_vert_shader("nebula");
		s.set_frag_shader("nebula");
		s.begin_shader();
		s.add_uniform_int("noise_tex", 0);
		s.add_uniform_float("noise_scale", 0.01);
		bind_3d_texture(get_noise_tex_3d(32, 4)); // RGBA noise
	}
	s.enable();
	enable_blend();
	glDepthMask(GL_FALSE); // no depth writing
}


void unebula::end_render(shader_t &s) {

	glDepthMask(GL_TRUE);
	disable_blend();
	s.disable();
}


void unebula::draw(point_d pos_, point const &camera, float max_dist, shader_t &s) const { // Note: new VFC here

	pos_ += pos;
	float const dist(p2p_dist(camera, pos_)), dist_scale(CLIP_TO_01(1.0f - 1.6f*(dist - radius)/max_dist));
	if (dist_scale <= 0.0) return; // too far away
	if (!univ_sphere_vis(pos_, radius)) return;
	float size_scale(1.0);
	move_in_front_of_far_clip(pos_, camera, size_scale, (dist + radius)); // distance to furthest point
	glPushMatrix();
	global_translate(pos_);
	uniform_scale(size_scale);
	colorRGBA mod_color[3];

	for (unsigned d = 0; d < 3; ++d) {
		mod_color[d] = color[d];
		mod_color[d].alpha *= ((draw_model == 1) ? 1.0 : 0.25*dist_scale);
	}
	if (s.is_setup()) { // assert setup?
		s.add_uniform_color("color1", mod_color[0]);
		s.add_uniform_color("color2", mod_color[1]);
		s.add_uniform_color("color3", mod_color[2]);
		s.add_uniform_vector3d("view_dir", (camera - pos_).get_norm()); // local object space
		s.add_uniform_float("radius", radius);
		s.add_uniform_float("offset", pos.x);
	}
	mod_color[0].do_glColor();
	points.front().set_state();
	glDrawArrays(GL_QUADS, 0, (unsigned)points.size());
	glPopMatrix();
}


