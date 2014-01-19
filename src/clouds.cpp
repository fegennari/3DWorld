// 3D World - Cloud/Nebula Generation and Drawing Code
// by Frank Gennari
// 3/10/12
#include "3DWorld.h"
#include "physics_objects.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "draw_utils.h"
#include "universe.h" // for unebula


bool     const USE_CLOUD_FBO    = 1;
unsigned const CLOUD_GEN_TEX_SZ = 1024;
unsigned const CLOUD_NUM_DIV = 32;


vector2d cloud_wind_pos(0.0, 0.0);
cloud_manager_t cloud_manager;

extern bool have_sun, no_sun_lpos_update;
extern int window_width, window_height, cloud_model, draw_model, display_mode, xoff, yoff, animate2, is_cloudy;
extern float CLOUD_CEILING, atmosphere, sun_rot, fticks, water_plane_z, zmin, zmax;
extern vector3d wind;
extern colorRGBA sun_color;


void draw_part_clouds(vector<particle_cloud> const &pc, colorRGBA const &color, bool zoomed);


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
	vector<cloud_t> clouds;

	if (calc_sun_light) {
		int last_src(0);

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
			light = max(0.5f, light);

			if (light_factor < 0.6) {
				float const blend(sqrt(5.0*(light_factor - 0.4)));
				light = light*blend + 0.25f*(1.0 - blend);
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
	draw_part_clouds(*this, WHITE, 1); // draw clouds
	camera_pos = orig_cpos;
	camera_pdu.valid = was_valid;
	set_red_only(0);

	if (!USE_CLOUD_FBO) { // render clouds to texture
		bind_2d_texture(cloud_tid);
		glReadBuffer(GL_BACK);
		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, xsize, ysize); // copy the frame buffer to the bound texture
	}

	// reset state
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	if (USE_CLOUD_FBO) disable_fbo();
	set_standard_viewport();
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
	//RESET_TIME;
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
		point const camera(get_camera_pos());
		cube_t const bcube(get_bcube());
		float const cloud_bot(bcube.d[2][0]), cloud_top(bcube.d[2][1]), cloud_xy(get_max_xy_extent());
		float const xy_exp((cloud_top - frustum_z)/(cloud_bot - frustum_z));
		//if (!camera_pdu.cube_visible(bcube)) return; // incorrect, and rarely returns

		create_texture(need_update);
		enable_flares(get_cloud_color(), 1); // texture will be overriden
		assert(cloud_tid);
		bind_2d_texture(cloud_tid);

		shader_t s;
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("cloud_billboard");
		s.begin_shader();
		s.add_uniform_int("tex0", 0);

		quad_batch_draw qbd;
		qbd.add_quad_dirs(point(camera.x, camera.y, cloud_top), vector3d(-xy_exp*cloud_xy, 0.0, 0.0), vector3d(0.0, xy_exp*cloud_xy, 0.0), WHITE, -plus_z);
		qbd.draw();
		s.end_shader();
		disable_flares();
	}
	else {
		draw_part_clouds(*this, get_cloud_color(), 1);
	}
	glEnable(GL_DEPTH_TEST);
	//PRINT_TIME("Clouds");
}


void draw_puffy_clouds(int order) {

	if (cloud_manager.is_inited() && (get_camera_pos().z > cloud_manager.get_z_plane()) != order) return;

	if (atmosphere < 0.01) {
		cloud_manager.clear();
	}
	else if (((display_mode & 0x40) != 0) ^ is_cloudy) { // key 7
		cloud_manager.draw();
	}
}


float get_cloud_zmax() {return get_camera_pos().z + max(zmax, CLOUD_CEILING);}


void set_cloud_uniforms(shader_t &s, unsigned tu_id) {

	select_multitex(NOISE_GEN_TEX, tu_id, 0);
	s.add_uniform_int("cloud_noise_tex", tu_id);
	s.add_uniform_vector2d("dxy", cloud_wind_pos);
}


void render_spherical_section(indexed_mesh_draw<vert_wrap_t> &imd, float size, float rval_inv, float z1, float z2) {

	point const camera(get_camera_pos());
	float const ndiv_inv(1.0/CLOUD_NUM_DIV);
	float xvals[CLOUD_NUM_DIV+1], yvals[CLOUD_NUM_DIV+1];
	imd.init(CLOUD_NUM_DIV, CLOUD_NUM_DIV);

	for (unsigned t = 0; t <= CLOUD_NUM_DIV; ++t) {
		float const angle(TWO_PI*(t*ndiv_inv));
		xvals[t] = cosf(angle);
		yvals[t] = sinf(angle);
	}
	for (unsigned r = 0; r <= CLOUD_NUM_DIV; ++r) { // spherical section
		float const radius(size*r*ndiv_inv), zval(z1 + (z2 - z1)*cos(PI_TWO*min(1.0f, radius*rval_inv)));

		for (unsigned t = 0; t <= CLOUD_NUM_DIV; ++t) {
			imd.set_vert(t, r, point((radius*xvals[t] + camera.x), (radius*yvals[t] + camera.y), zval));
		}
	}
	imd.render();
}


// not a plane, but a spherical section
void draw_cloud_plane(float terrain_zmin, bool reflection_pass) {

	float const cloud_rel_vel = 1.0; // relative cloud velocity compared to camera velocity (0: clouds follow the camera, 1: clouds are stationary)
	float const size(camera_pdu.far_), rval(0.94*size), rval_inv(1.0/rval); // extends to at least the far clipping plane
	float const z1(zmin), z2(get_cloud_zmax()), ndiv_inv(1.0/CLOUD_NUM_DIV);
	point const camera(get_camera_pos()), world_pos(camera + vector3d((xoff2-xoff)*DX_VAL, (yoff2-yoff)*DY_VAL, 0.0));
	vector3d const offset(-camera + cloud_rel_vel*world_pos);
	colorRGBA const cloud_color(get_cloud_color());
	static indexed_mesh_draw<vert_wrap_t> imd;
	shader_t s;
	glDepthMask(GL_FALSE);

	if (animate2) {
		cloud_wind_pos.x -= fticks*wind.x;
		cloud_wind_pos.y -= fticks*wind.y;
	}

	// draw a static textured upper cloud layer
	if (1) {
		s.set_vert_shader("texture_gen.part+no_lighting_texture_gen");
		s.set_frag_shader("simple_texture");
		s.begin_shader();
		s.add_uniform_int("tex0", 0);
		enable_blend();
		select_texture(CLOUD_RAW_TEX, 0);
		float const xy_scale(0.02), z_offset(0.01*size);
		float const sx(1.0/(1.0 + min(2.0, 0.5*fabs(wind.x))));
		float const sy(1.0/(1.0 + min(2.0, 0.5*fabs(wind.y))));
		vector3d const offset2(xy_scale*(-camera + 0.7*cloud_rel_vel*world_pos));
		setup_texgen(sx*xy_scale, sy*xy_scale, sx*offset2.x, sy*offset2.y, 0.0);
		colorRGBA cloud_color2(cloud_color);
		cloud_color2.alpha *= (is_cloudy ? 1.0 : 0.5);
		cloud_color2.do_glColor();
		render_spherical_section(imd, size, rval_inv, z1+z_offset, z2+z_offset);
		disable_textures_texgen();
		disable_blend();
		s.end_shader();
	}

	// draw a plane at terrain_zmin to properly blend the fog
	if (!reflection_pass && glIsEnabled(GL_FOG)) {
		colorRGBA fog_color;
		glGetFloatv(GL_FOG_COLOR, (float *)&fog_color);
		fog_color.do_glColor();
		s.begin_color_only_shader();
		imd.render_z_plane(-size, -size, size, size, (terrain_zmin - SMALL_NUMBER), CLOUD_NUM_DIV, CLOUD_NUM_DIV);
		s.end_shader();
	}

	// draw clouds
	if ((display_mode & 0x40) == 0) { // on by default
		setup_tt_fog_pre(s);
		s.set_prefix("#define NUM_OCTAVES 8", 1); // FS
		s.set_bool_prefix("underwater_atten", (glIsEnabled(GL_FOG) != 0), 1); // FS
		s.set_vert_shader("water_fog.part+clouds");
		s.set_frag_shader("linear_fog.part+perlin_clouds.part*+clouds");
		s.begin_shader();
		setup_tt_fog_post(s);
		s.add_uniform_float("water_plane_z", zmin);
		s.add_uniform_float("cloud_scale", (is_cloudy ? 1.0 : 0.5));
		s.add_uniform_vector3d("camera_pos", camera);
		s.add_uniform_vector3d("sun_pos", get_sun_pos());
		s.add_uniform_color("sun_color", sun_color);
		set_cloud_uniforms(s, 0);
		s.add_uniform_vector3d("cloud_offset", offset);
		enable_blend();
		cloud_color.do_glColor();
		render_spherical_section(imd, size, rval_inv, z1, z2);
		disable_blend();
		s.end_shader();
	}
	glDepthMask(GL_TRUE);
}


// *** nebula code ***


void move_in_front_of_far_clip(point_d &pos, point const &camera, float &size, float dist, float dscale);


/*static*/ colorRGBA volume_part_cloud::gen_color(rand_gen_t &rgen) {

	return colorRGBA(rgen.rand_uniform(0.3, 1.0), rgen.rand_uniform(0.1, 0.5), rgen.rand_uniform(0.2, 0.9), 1.0);
}


void volume_part_cloud::gen_pts(float radius) {

	unsigned ix(0);
	points.resize(4*13);

	for (int z = -1; z <= 1; ++z) {
		for (int y = -1; y <= 1; ++y) {
			for (int x = -1; x <= 1; ++x) {
				if (x == 0 && y == 0 && z == 0) continue; // zero normal
				int v(1); if (x) {v *= x;} if (y) {v *= y;} if (z) {v *= z;}
				if (v < 0) continue; // skip mirrored normals
				vector3d const normal(vector3d(x, y, z).get_norm());
				vector3d vab[2];
				get_ortho_vectors(normal, vab);

				for (unsigned j = 0; j < 4; ++j) { // Note: quads will extend beyond radius, but will be rendered as alpha=0 outside radius
					points[ix+j].v = radius*(((j>>1) ? 1.0 : -1.0)*vab[0] + (((j&1)^(j>>1)) ? 1.0 : -1.0)*vab[1]);
					points[ix+j].set_norm(normal);
				}
				ix += 4;
			}
		}
	}
	assert(ix == points.size());
}


/*static*/ void volume_part_cloud::shader_setup(shader_t &s, unsigned noise_ncomp) {

	assert(noise_ncomp == 1 || noise_ncomp == 4);
	bind_3d_texture(get_noise_tex_3d(32, noise_ncomp));
	if (s.is_setup()) return; // nothing else to do
	s.set_prefix("#define NUM_OCTAVES 5", 1); // FS
	s.set_int_prefix("noise_ncomp", noise_ncomp, 1); // FS
	s.set_bool_prefix("line_mode", (draw_model == 1), 1); // FS
	s.set_vert_shader("nebula");
	s.set_frag_shader("nebula");
	s.begin_shader();
	s.add_uniform_int("noise_tex", 0);
}


void volume_part_cloud::draw_quads() const {

	glDepthMask(GL_FALSE); // no depth writing
	draw_verts(points, GL_QUADS);
	glDepthMask(GL_TRUE);
}


void unebula::gen(float range, ellipsoid_t const &bounds) {

	// Note: bounds is not currently used, but it can be used to scale the nebula to the galaxy's ellipsoid (but requires some transforms in the shader)
	rand_gen_t rgen;
	rgen.set_state(rand2(), rand2());
	radius = rgen.rand_uniform(0.1, 0.15)*range;
	UNROLL_3X(color[i_] = gen_color(rgen);)
	gen_pts(radius);
}


void unebula::draw(point_d pos_, point const &camera, float max_dist, shader_t &s) const { // Note: new VFC here

	pos_ += pos;
	float const dist(p2p_dist(camera, pos_)), dist_scale(CLIP_TO_01(1.0f - 1.6f*(dist - radius)/max_dist));
	if (dist_scale <= 0.0) return; // too far away
	if (!univ_sphere_vis(pos_, radius)) return;
	float size_scale(1.0);
	move_in_front_of_far_clip(pos_, camera, size_scale, (dist + radius), 1.5); // distance to furthest point
	glPushMatrix();
	global_translate(pos_);
	uniform_scale(size_scale);
	colorRGBA mod_color[3];

	for (unsigned d = 0; d < 3; ++d) {
		mod_color[d] = color[d];
		mod_color[d].alpha *= ((draw_model == 1) ? 1.0 : 0.3*dist_scale);
	}
	shader_setup(s, 4); // RGBA noise
	s.enable();
	s.add_uniform_float("noise_scale", 0.01);
	s.add_uniform_color("color1i", mod_color[0]);
	s.add_uniform_color("color1o", mod_color[0]);
	s.add_uniform_color("color2i", mod_color[1]);
	s.add_uniform_color("color2o", mod_color[1]);
	s.add_uniform_color("color3i", mod_color[2]);
	s.add_uniform_color("color3o", mod_color[2]);
	s.add_uniform_float("radius",  radius);
	s.add_uniform_float("offset",  pos.x); // used as a hash
	s.add_uniform_vector3d("view_dir", (camera - pos_).get_norm()); // local object space
	mod_color[0].do_glColor();
	enable_blend();
	draw_quads();
	disable_blend();
	s.disable();
	glPopMatrix();
}


