// 3D World - Cloud Generation and Drawing Code
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "physics_objects.h"
#include "shaders.h"
#include "gl_ext_arb.h"


unsigned const CLOUD_GEN_TEX_SZ = 1024;


cloud_manager_t cloud_manager;

extern bool have_sun, no_sun_lpos_update;
extern int window_width, window_height, cloud_model;
extern float CLOUD_CEILING, atmosphere, sun_rot;


void draw_part_cloud(vector<particle_cloud> const &pc, colorRGBA const color, bool zoomed);


struct cloud_t {
	unsigned b, e;
	point p;
	float r;
	cloud_t(unsigned b_=0, unsigned e_=0, point const &p_=all_zeros, float r_=0.0)
		: b(b_), e(e_), p(p_), r(r_) {}
};


void cloud_manager_t::create_clouds() { // 3D cloud puffs

	unsigned const NCLOUDS = 10;
	unsigned const NPARTS  = 1000;
	clear();
	srand(123);

	for (unsigned c = 0; c < NCLOUDS; ++c) {
		point const center(4.0*X_SCENE_SIZE*signed_rand_float(),
			               4.0*Y_SCENE_SIZE*signed_rand_float(),
						   (ztop + CLOUD_CEILING + Z_SCENE_SIZE*rand_uniform(0.25, 0.75)));
		point const bounds(X_SCENE_SIZE*rand_uniform(1.0, 2.0),
			               Y_SCENE_SIZE*rand_uniform(1.0, 2.0),
						   Z_SCENE_SIZE*rand_uniform(0.4, 0.8));
		unsigned const nparts(rand()%(NPARTS/2) + NPARTS/2);
		unsigned const ix(size());
		resize(ix + nparts);

		for (unsigned p = 0; p < nparts; ++p) {
			point pos(signed_rand_vector_spherical(1.0));

			for (unsigned i = 0; i < 3; ++i) {
				pos[i] *= bounds[i];
			}
			if (pos.z < 0.0) pos.z *= 0.5; // compressed on the bottom
			pos += center;
			float const radius(0.045*(X_SCENE_SIZE + Y_SCENE_SIZE)*rand_uniform(0.5, 1.0));
			float const density(rand_uniform(0.05, 0.12));
			(*this)[ix + p].gen(pos, WHITE, zero_vector, radius, density, 0.0, 0.0, -((int)c+2), 0, 0, 1, 1); // no lighting
		}
	}
}


void cloud_manager_t::update_lighting() {

	RESET_TIME;
	point const sun_pos(get_sun_pos());
	bool const calc_sun_light(have_sun && light_factor > 0.4);
	unsigned const num_clouds(size());
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


bool cloud_manager_t::create_texture(bool force_recreate) {

	RESET_TIME;
	unsigned const xsize(min(CLOUD_GEN_TEX_SZ, (unsigned)window_width)), ysize(min(CLOUD_GEN_TEX_SZ, (unsigned)window_height));

	if (txsize != xsize || tysize != ysize) {
		free_textures();
		txsize = xsize;
		tysize = ysize;
	}
	if (cloud_tid && !force_recreate) return 0; // nothing to do
	
	if (!cloud_tid) {
		setup_texture(cloud_tid, GL_MODULATE, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, xsize, ysize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	}
	assert(glIsTexture(cloud_tid));
	check_gl_error(800);

	//enable_fbo(fbo_id, cloud_tid, 0);

	glViewport(0, 0, xsize, ysize);
	glClearColor(1.0, 1.0, 1.0, 1.0); // white
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	// setup projection matrix
	cube_t const bcube(get_bcube());
	float const dx(max(fabs(bcube.d[0][0]), fabs(bcube.d[0][1]))), dy(max(fabs(bcube.d[1][0]), fabs(bcube.d[1][1])));
	float const cloud_bot(bcube.d[2][0]), cloud_top(bcube.d[2][1]);
	float const cloud_xy(max(dx, dy)), scene_xy(max(X_SCENE_SIZE, Y_SCENE_SIZE)), angle(atan2(cloud_xy, cloud_bot));
	float const z1(min(zbottom, czmin)), frustum_z(z1 - scene_xy*(cloud_bot - z1)/(cloud_xy - scene_xy));
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
	bool const was_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable view frustum culling
	draw_part_cloud(*this, WHITE, 1); // draw clouds
	camera_pdu.valid = was_valid;
	set_red_only(0);

	// render clouds to texture
	glBindTexture(GL_TEXTURE_2D, cloud_tid);
	glReadBuffer(GL_BACK);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, xsize, ysize); // copy the frame buffer to the bound texture

	// reset state
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glViewport(0, 0, window_width, window_height);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	check_gl_error(801);
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
	if (empty()) create_clouds();
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
	if (cloud_model == 1) {
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
		cube_t const bcube(get_bcube());
		
		for (unsigned d = 0; d < 2; ++d) { // render the bottom face of bcube
			for (unsigned e = 0; e < 2; ++e) {
				glTexCoord2f(float(d^e^1), float(d));
				point(bcube.d[0][d^e], bcube.d[1][d], bcube.d[2][0]).do_glVertex();
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


