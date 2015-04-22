// 3D World - Cloud/Nebula Generation and Drawing Code
// by Frank Gennari
// 3/10/12
#include "3DWorld.h"
#include "physics_objects.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "draw_utils.h"
#include "universe.h" // for unebula
#include "cobj_bsp_tree.h"


unsigned const CLOUD_GEN_TEX_SZ = 1024;
unsigned const CLOUD_NUM_DIV = 32;


vector2d cloud_wind_pos(0.0, 0.0);
cloud_manager_t cloud_manager;

extern bool have_sun, no_sun_lpos_update, fog_enabled;
extern int window_width, window_height, cloud_model, draw_model, display_mode, xoff, yoff, animate2, is_cloudy;
extern float CLOUD_CEILING, NEAR_CLIP, FAR_CLIP, atmosphere, sun_rot, fticks, cloud_cover, zmin, zmax;
extern vector3d wind;
extern colorRGBA sun_color, cur_fog_color;


void draw_part_clouds(vector<particle_cloud> const &pc, int tid);


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


class cloud_bvh_t : public cobj_tree_sphere_t {

	cloud_manager_t &mgr;

public:
	cloud_bvh_t(cloud_manager_t &mgr_) : mgr(mgr_) {}

	void setup(bool verbose) {
		objects.resize(mgr.size());
		for (unsigned i = 0; i < mgr.size(); ++i) {objects[i] = sphere_with_id_t(mgr[i].pos, mgr[i].radius, i);}
		build_tree_top(verbose);
	}
	float calc_light_value(point const &pos, point const &sun_pos) const {
		vector3d const v1(sun_pos - pos);
		float const dist_sq(v1.mag_sq());
		vector3d const v1n(v1/dist_sq);
		float light(1.0); // start off fully lit
		node_ix_mgr nixm(nodes, pos, sun_pos);
		unsigned const num_nodes((unsigned)nodes.size());

		for (unsigned nix = 0; nix < num_nodes;) {
			tree_node const &n(nodes[nix]);
			if (!nixm.check_node(nix)) continue;

			for (unsigned i = n.start; i < n.end; ++i) { // check leaves
				particle_cloud const &pc(mgr[objects[i].id]);
				vector3d const v2(sun_pos, pc.pos);
				if (v2.mag_sq() > dist_sq) continue; // further from the sun
				float const dotp(dot_product(v1, v2));
				float const dsq((dotp > dist_sq) ? p2p_dist_sq(v1, v2) : (v2 - v1n*dotp).mag_sq());
				if (dsq > pc.radius*pc.radius) continue; // no intersection
				float const alpha(2.0*pc.base_color.alpha*pc.density*((pc.radius - sqrt(dsq))/pc.radius));
				light *= 1.0 - CLIP_TO_01(alpha);
			} // for i
		} // for nix
		return light;
	}
};


void cloud_manager_t::update_lighting() {

	RESET_TIME;
	point const sun_pos(get_sun_pos());
	bool const calc_sun_light(have_sun && light_factor > 0.4);
	unsigned const num_clouds((unsigned)size());

	if (!calc_sun_light) {
		for (unsigned i = 0; i < num_clouds; ++i) {
			particle_cloud &pc((*this)[i]);
			pc.darkness   = 0.5; // night time sky
			pc.base_color = WHITE;
			apply_red_sky(pc.base_color);
		}
		return;
	}
	cloud_bvh_t bvh(*this);
	bvh.setup(0);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)num_clouds; ++i) {
		particle_cloud &pc((*this)[i]);
		float light(max(0.5f, bvh.calc_light_value(pc.pos, sun_pos)));

		if (light_factor < 0.6) {
			float const blend(sqrt(5.0*(light_factor - 0.4)));
			light = light*blend + 0.25f*(1.0 - blend);
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
	unsigned const xsize(CLOUD_GEN_TEX_SZ), ysize(CLOUD_GEN_TEX_SZ);

	if (txsize != xsize || tysize != ysize) {
		free_textures();
		txsize = xsize;
		tysize = ysize;
	}
	if (cloud_tid && !force_recreate) return 0; // nothing to do
	
	if (!cloud_tid) {
		setup_texture(cloud_tid, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, xsize, ysize, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	}
	assert(glIsTexture(cloud_tid));
	check_gl_error(800);
	enable_fbo(fbo_id, cloud_tid, 0);
	check_gl_error(801);

	glViewport(0, 0, xsize, ysize);
	glClearColor(1.0, 1.0, 1.0, 1.0); // white
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	fgMatrixMode(FG_PROJECTION);
	fgPushMatrix();
	fgLoadIdentity();

	// setup projection matrix
	cube_t const bcube(get_bcube());
	float const cloud_bot(bcube.d[2][0]), cloud_top(bcube.d[2][1]), cloud_xy(get_max_xy_extent());
	float const scene_xy(max(X_SCENE_SIZE, Y_SCENE_SIZE)), angle(atan2(cloud_xy, cloud_bot)), z1(min(zbottom, czmin));
	frustum_z = z1 - scene_xy*(cloud_bot - z1)/(cloud_xy - scene_xy);
	//pos_dir_up const pdu(get_pt_cube_frustum_pdu(get_camera_pos(), bcube, 1));
	//pos_dir_up const pdu(all_zeros, plus_z, plus_x, tanf(angle)*SQRT2, sinf(angle), NEAR_CLIP, FAR_CLIP, 1.0);
	//fgPerspective(2.0*angle/TO_RADIANS, 1.0, cloud_bot-frustum_z, cloud_top+(cloud_top - cloud_bot)-frustum_z);
	fgPerspective(2.0*angle/TO_RADIANS, 1.0, NEAR_CLIP, FAR_CLIP);
	fgMatrixMode(FG_MODELVIEW);
	fgPushMatrix();
	fgLoadIdentity();
	vector3d const up_dir(plus_y);
	point const origin(0.0, 0.0, frustum_z), center(0.0, 0.0, cloud_bot);
	fgLookAt(origin.x, origin.y, origin.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);

	set_red_only(1);
	point const orig_cpos(camera_pos);
	bool const was_valid(camera_pdu.valid);
	camera_pdu.valid = 0; // disable view frustum culling
	camera_pos = origin;
	shader_t s;
	s.begin_simple_textured_shader(0.01);
	draw_part_clouds(*this, SMOKE_PUFF_TEX); // draw clouds
	s.end_shader();
	camera_pos = orig_cpos;
	camera_pdu.valid = was_valid;
	set_red_only(0);
	// reset state
	fgMatrixMode(FG_PROJECTION);
	fgPopMatrix();
	fgMatrixMode(FG_MODELVIEW);
	fgPopMatrix();
	disable_fbo();
	set_standard_viewport();
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

	// WRITE: wind moves clouds

	// light source code
	static bool had_sun(0);
	static float last_sun_rot(0.0);
	bool const need_update(!no_sun_lpos_update && (sun_rot != last_sun_rot || have_sun != had_sun));
	int const tid(SMOKE_PUFF_TEX);
	set_multisample(0);
	glDisable(GL_DEPTH_TEST);
	shader_t s;

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
		enable_flares(tid); // texture will be overriden
		assert(cloud_tid);
		bind_2d_texture(cloud_tid);

		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("cloud_billboard");
		s.begin_shader();
		s.add_uniform_int("tex0", 0);

		quad_batch_draw qbd;
		qbd.add_quad_dirs(point(camera.x, camera.y, cloud_top), vector3d(-xy_exp*cloud_xy, 0.0, 0.0), vector3d(0.0, xy_exp*cloud_xy, 0.0), get_cloud_color(), -plus_z);
		qbd.draw();
		disable_flares();
	}
	else {
		s.begin_simple_textured_shader(0.01);
		draw_part_clouds(*this, tid);
	}
	s.end_shader();
	glEnable(GL_DEPTH_TEST);
	set_multisample(1);
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

	select_multitex(NOISE_GEN_TEX, tu_id);
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


indexed_mesh_draw<vert_wrap_t> cloud_imd;

void free_cloud_context() {cloud_imd.free_context();}


vector3d get_cloud_offset(float rel_vel_scale) {

	float const cloud_rel_vel = 1.0; // relative cloud velocity compared to camera velocity (0: clouds follow the camera, 1: clouds are stationary)
	point const camera(get_camera_pos()), world_pos(camera + vector3d((xoff2-xoff)*DX_VAL, (yoff2-yoff)*DY_VAL, 0.0));
	return -camera + rel_vel_scale*cloud_rel_vel*world_pos;
}


// not a plane, but a spherical section
void draw_cloud_planes(float terrain_zmin, bool reflection_pass, bool draw_ceil, bool draw_floor) {

	shader_t s;
	float const size(camera_pdu.far_);

	// draw a plane at terrain_zmin to properly blend the fog (needs to be first when camera is above the clouds)
	if (draw_floor && !reflection_pass && fog_enabled) {
		glDepthMask(GL_FALSE);
		s.begin_color_only_shader(cur_fog_color);
		cloud_imd.render_z_plane(-size, -size, size, size, (terrain_zmin - SMALL_NUMBER), CLOUD_NUM_DIV, CLOUD_NUM_DIV);
		s.end_shader();
		glDepthMask(GL_TRUE);
	}
	if (!draw_ceil) return;

	// track the min values of terrain_zmin seen and reset when zmin changes (mesh change)
	static float prev_zmin(0.0), min_terrain_zmin(0.0);
	if (zmin != prev_zmin) {prev_zmin = zmin; min_terrain_zmin = terrain_zmin;}
	else {min_terrain_zmin = min(min_terrain_zmin, terrain_zmin);}

	float const cloud_rel_vel = 1.0; // relative cloud velocity compared to camera velocity (0: clouds follow the camera, 1: clouds are stationary)
	float const rval(0.94*size), rval_inv(1.0/rval); // extends to at least the far clipping plane
	float const cloud_z(get_tt_cloud_level()); // halfway between the top of the mountains and the end of the atmosphere
	float const z1(min(zmin, min_terrain_zmin)), z2(min(cloud_z, get_cloud_zmax())), ndiv_inv(1.0/CLOUD_NUM_DIV);
	vector3d const offset(get_cloud_offset(1.0));
	colorRGBA const cloud_color(get_cloud_color());

	if (animate2) {
		cloud_wind_pos.x -= fticks*wind.x;
		cloud_wind_pos.y -= fticks*wind.y;
	}

	// draw a static textured upper cloud layer
	glDepthMask(GL_FALSE);

	if (1) {
		s.set_vert_shader("texture_gen.part+no_lighting_texture_gen");
		s.set_frag_shader("simple_texture");
		s.begin_shader();
		s.add_uniform_int("tex0", 0);
		enable_blend();
		select_texture(CLOUD_RAW_TEX);
		float const xy_scale(0.02), z_offset(0.01*size);
		float const sx(1.0/(1.0 + min(2.0, 0.5*fabs(wind.x))));
		float const sy(1.0/(1.0 + min(2.0, 0.5*fabs(wind.y))));
		vector3d const offset2(xy_scale*get_cloud_offset(0.7));
		setup_texgen(sx*xy_scale, sy*xy_scale, sx*offset2.x, sy*offset2.y, 0.0, s, 0);
		colorRGBA cloud_color2(cloud_color);
		cloud_color2.alpha *= (is_cloudy ? 1.0 : 0.5);
		s.set_cur_color(cloud_color2);
		render_spherical_section(cloud_imd, size, rval_inv, z1+z_offset, z2+z_offset);
		disable_blend();
		s.end_shader();
	}

	// draw clouds
	if ((display_mode & 0x40) == 0) { // on by default
		setup_tt_fog_pre(s);
		s.set_prefix("#define NUM_OCTAVES 8", 1); // FS
		s.set_bool_prefix("underwater_atten", fog_enabled, 1); // FS
		s.set_vert_shader("water_fog.part*+clouds");
		s.set_frag_shader("linear_fog.part+perlin_clouds.part*+clouds"); // Note: could also apply water fog in fragment shader
		s.begin_shader();
		setup_tt_fog_post(s);
		s.add_uniform_float("water_plane_z", z1);
		s.add_uniform_float("cloud_plane_z", z2);
		s.add_uniform_float("cloud_scale", (is_cloudy ? 1.0 : (cloud_cover + 0.5*(1.0 - cloud_cover))));
		s.add_uniform_vector3d("camera_pos", get_camera_pos());
		s.add_uniform_vector3d("sun_pos", get_sun_pos()); // no sun at night?
		s.add_uniform_color("sun_color", sun_color);
		set_cloud_uniforms(s, 0);
		s.add_uniform_vector3d("cloud_offset", offset);
		enable_blend();
		s.set_cur_color(cloud_color);
		render_spherical_section(cloud_imd, size, rval_inv, z1, z2);
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


/*static*/ vector<volume_part_cloud::vert_type_t> volume_part_cloud::unscaled_points;

/*static*/ void volume_part_cloud::cacl_unscaled_points() {

	unsigned ix(0);
	unscaled_points.resize(4*13);

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
					unscaled_points[ix+j].v = ((j>>1) ? 1.0 : -1.0)*vab[0] + (((j&1)^(j>>1)) ? 1.0 : -1.0)*vab[1];
					unscaled_points[ix+j].set_norm(normal);
				}
				ix += 4;
			}
		}
	}
	assert(ix == unscaled_points.size());
}


void volume_part_cloud::gen_pts(float radius) {

	if (unscaled_points.empty()) {cacl_unscaled_points();}
	points = unscaled_points;
	for (unsigned i = 0; i < points.size(); ++i) {points[i].v *= radius;}
}


void vpc_shader_t::cache_locs() {

	ns_loc  = get_uniform_loc("noise_scale");
	c1i_loc = get_uniform_loc("color1i");
	c1o_loc = get_uniform_loc("color1o");
	c2i_loc = get_uniform_loc("color2i");
	c2o_loc = get_uniform_loc("color2o");
	c3i_loc = get_uniform_loc("color3i");
	c3o_loc = get_uniform_loc("color3o");
	rad_loc = get_uniform_loc("radius");
	off_loc = get_uniform_loc("offset");
	vd_loc  = get_uniform_loc("view_dir");
}


/*static*/ void volume_part_cloud::shader_setup(vpc_shader_t &s, unsigned noise_ncomp) {

	assert(noise_ncomp == 1 || noise_ncomp == 4);
	bind_3d_texture(get_noise_tex_3d(32, noise_ncomp));
	if (s.is_setup()) return; // nothing else to do
	s.set_prefix("#define NUM_OCTAVES 5", 1); // FS
	s.set_prefix("#define RIDGED_NOISE",  1); // FS
	s.set_int_prefix("noise_ncomp", noise_ncomp, 1); // FS
	s.set_bool_prefix("line_mode", (draw_model == 1), 1); // FS
	s.set_vert_shader("nebula");
	s.set_frag_shader("nebula");
	s.begin_shader();
	s.add_uniform_int("noise_tex", 0);
	s.cache_locs();
}


void volume_part_cloud::draw_quads() const {

	glDepthMask(GL_FALSE); // no depth writing
	draw_quad_verts_as_tris(points);
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


void unebula::draw(point_d pos_, point const &camera, float max_dist, vpc_shader_t &s) const { // Note: new VFC here

	pos_ += pos;
	float const dist(p2p_dist(camera, pos_)), dist_scale(CLIP_TO_01(1.0f - 1.6f*(dist - radius)/max_dist));
	if (dist_scale <= 0.0) return; // too far away
	if (!univ_sphere_vis(pos_, radius)) return;
	float size_scale(1.0);
	move_in_front_of_far_clip(pos_, camera, size_scale, (dist + radius), 1.5); // distance to furthest point
	fgPushMatrix();
	global_translate(pos_);
	uniform_scale(size_scale);
	colorRGBA mod_color[3];

	for (unsigned d = 0; d < 3; ++d) {
		mod_color[d] = color[d];
		mod_color[d].alpha *= ((draw_model == 1) ? 1.0 : 0.3*dist_scale);
	}
	shader_setup(s, 4); // RGBA noise
	s.enable();
	s.set_uniform_float(s.ns_loc, 0.01);
	s.set_uniform_color(s.c1i_loc, mod_color[0]);
	s.set_uniform_color(s.c1o_loc, mod_color[0]);
	s.set_uniform_color(s.c2i_loc, mod_color[1]);
	s.set_uniform_color(s.c2o_loc, mod_color[1]);
	s.set_uniform_color(s.c3i_loc, mod_color[2]);
	s.set_uniform_color(s.c3o_loc, mod_color[2]);
	s.set_uniform_float(s.rad_loc, radius);
	s.set_uniform_float(s.off_loc, pos.x); // used as a hash
	s.set_uniform_vector3d(s.vd_loc, (camera - pos_).get_norm()); // local object space
	s.set_cur_color(mod_color[0]);
	enable_blend();
	draw_quads();
	disable_blend();
	s.disable();
	fgPopMatrix();
}


