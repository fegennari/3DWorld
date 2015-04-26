// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "3DWorld.h"
#include "collision_detect.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "model3d.h"
#include "shadow_map.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>


bool scene_smap_vbo_invalid(0), voxel_shadows_updated(0);
unsigned shadow_map_sz(0);
pos_dir_up orig_camera_pdu;

extern bool snow_shadows;
extern int window_width, window_height, animate2, display_mode, tree_mode, ground_effects_level, num_trees, camera_coll_id;
extern unsigned enabled_lights;
extern float NEAR_CLIP, tree_deadness, vegetation;
extern vector<shadow_sphere> shadow_objs;
extern coll_obj_group coll_objects;
extern platform_cont platforms;

void draw_trees(bool shadow_only=0);


cube_t get_scene_bounds() {
	return cube_t(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, min(zbottom, czmin), max(ztop, czmax));
}


struct ground_mode_smap_data_t : public smap_data_t {

	bool last_has_dynamic;

	ground_mode_smap_data_t(unsigned tu_id_) : smap_data_t(tu_id_), last_has_dynamic(0) {}
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
};

ground_mode_smap_data_t smap_data[NUM_LIGHT_SRC] = {ground_mode_smap_data_t(6), ground_mode_smap_data_t(7)};


class smap_vertex_cache_t : public vbo_wrap_t {

	unsigned num_verts1, num_verts2;

public:
	vector<vert_wrap_t> dverts;

	smap_vertex_cache_t() : num_verts1(0), num_verts2(0) {}
	void end_block1(unsigned size) {num_verts1 = size;}

	void upload(vector<vert_wrap_t> const &verts) {
		if (!verts.empty()) {create_and_upload(verts);}
		assert(vbo_valid());
		num_verts2 = verts.size();
	}
	void render() const {
		if (num_verts2 == 0) return; // empty
		shader_t s;
		s.begin_color_only_shader(WHITE); // Note: color is likely unused
		assert(num_verts1 <= num_verts2);
		pre_render();
		if (num_verts1 > 0) {draw_verts<vert_wrap_t>(NULL, num_verts1, GL_TRIANGLES);}

		if (num_verts2 > num_verts1) {
			glEnable(GL_CULL_FACE);
			//glCullFace(GL_FRONT); // faster, but artifacts at surface intersections such as stairs, especially with more z bias
			draw_verts<vert_wrap_t>(NULL, (num_verts2 - num_verts1), GL_TRIANGLES, num_verts1);
			//glCullFace(GL_BACK);
			glDisable(GL_CULL_FACE);
		}
		post_render();
		s.end_shader();
	}
	void render_dynamic() {
		draw_and_clear_verts(dverts, GL_TRIANGLES);
	}
	void free() {
		clear_vbo();
		num_verts1 = num_verts2 = 0;
		dverts.clear();
	}
};

smap_vertex_cache_t smap_vertex_cache;


bool shadow_map_enabled() {return (shadow_map_sz > 0);}

float approx_pixel_width() {
	return 0.5*sqrt(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE) / shadow_map_sz;
}

int get_smap_ndiv(float radius) {
	// dynamic based on distance(camera, line(lpos, scene_center))?
	return min(N_SPHERE_DIV, max(4, int(0.5*radius/approx_pixel_width())));
}

void free_smap_vbo() {smap_vertex_cache.free();}


xform_matrix get_texture_matrix(xform_matrix const &camera_mv_matrix) {
	
	// This matrix transforms every coordinate {x,y,z} to {x,y,z}* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	glm::mat4 const bias(	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
		0.5, 0.5, 0.5, 1.0);
	return (bias * fgGetPJM() * fgGetMVM() * glm::affineInverse((glm::mat4)camera_mv_matrix));
}


void smap_data_state_t::free_gl_state() {
	free_texture(tid);
	free_fbo(fbo_id);
}

bool smap_data_t::set_smap_shader_for_light(shader_t &s, int light, xform_matrix const *const mvm) const {

	if (!shadow_map_enabled() || !is_light_enabled(light)) return 0;
	point lpos; // unused
	bool const light_valid(light_valid(0xFF, light, lpos));
	s.add_uniform_int  (append_ix(string("sm_tex"),   light, 0), tu_id);
	s.add_uniform_float(append_ix(string("sm_scale"), light, 0), (light_valid ? 1.0 : 0.0));
	xform_matrix tm(texture_matrix);
	if (mvm) {tm *= (*mvm) * glm::affineInverse((glm::mat4)fgGetMVM());} // Note: works for translate, but not scale?
	s.add_uniform_matrix_4x4(append_ix(string("smap_matrix"), light, 0), tm.get_ptr(), 0);
	set_active_texture(tu_id);

	if (light_valid) { // otherwise, we know that sm_scale will be 0.0 and we won't do the lookup
		assert(tid > 0);
		bind_2d_texture(tid);
	}
	else {
		select_texture(WHITE_TEX); // default white texture
	}
	set_active_texture(0);
	return 1;
}


void set_smap_shader_for_all_lights(shader_t &s, float z_bias) {

	s.add_uniform_float("z_bias", z_bias);

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		smap_data[l].set_smap_shader_for_light(s, l);
	}
}


pos_dir_up get_pt_cube_frustum_pdu(point const &pos, cube_t const &bounds, bool set_matrix) {

	point const center(bounds.get_cube_center());
	vector3d const light_dir((center - pos).get_norm()); // almost equal to lpos (point light)
	float const dist(p2p_dist(pos, center));
	vector3d up_dir(zero_vector);
	up_dir[get_min_dim(light_dir)] = 1.0;
	point corners[8];
	get_cube_corners(bounds.d, corners);
	float const radius(bounds.get_bsphere_radius());

	// tighter bounds / higher quality / slower
	vector3d dirs[2]; // x, y (up)
	orthogonalize_dir(up_dir, light_dir, dirs[1], 1);
	cross_product(light_dir, dirs[1], dirs[0]);
	float rx(0.0), ry(0.0);

	for (unsigned i = 0; i < 8; ++i) {
		vector3d const delta((corners[i] - pos).get_norm());
		rx = max(rx, fabs(dot_product(dirs[0], delta)));
		ry = max(ry, fabs(dot_product(dirs[1], delta)));
	}
	float const frustum_skew_val(1.0 + 0.5*(bounds.d[2][1] - bounds.d[2][0])/dist);
	float const angle(atan2(frustum_skew_val*ry, 1.0f)), aspect(rx/ry);
	pos_dir_up const pdu(pos, light_dir, up_dir, angle, max(NEAR_CLIP, dist-radius), dist+radius, aspect);

	if (set_matrix) {
		fgMatrixMode(FG_PROJECTION);
		fgLoadIdentity();
		fgPerspective(2.0*angle/TO_RADIANS, aspect, pdu.near_, pdu.far_);
		fgMatrixMode(FG_MODELVIEW);
		fgLoadIdentity();
		fgLookAt(pos.x, pos.y, pos.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);
	}
	return pdu;
}


void draw_scene_bounds_and_light_frustum(point const &lpos) {

	// draw scene bounds
	shader_t s;
	enable_blend();
	s.begin_color_only_shader(colorRGBA(1.0, 1.0, 1.0, 0.25)); // white
	draw_simple_cube(get_scene_bounds(), 0);

	// draw light frustum
	s.begin_color_only_shader(colorRGBA(1.0, 1.0, 0.0, 0.25)); // yellow
	get_pt_cube_frustum_pdu(lpos, get_scene_bounds(), 0).draw_frustum();
	disable_blend();
	s.end_shader();
}


void set_shadow_tex_params() {

	// This is to allow usage of textureProj function in the shader
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
}


bool no_sparse_smap_update() {

	if (world_mode != WMODE_GROUND) return 0;
	bool const leaf_wind(num_trees > 0 && (display_mode & 0x0100) != 0 && (tree_mode & 1) != 0 && tree_deadness < 1.0 && vegetation > 0.0);
	if (leaf_wind || !shadow_objs.empty() || !platforms.empty()) return 1;
	//return !coll_objects.drawn_ids.empty();
	//return !coll_objects.dynamic_ids.empty();
	return 0;
}


bool smap_data_t::needs_update(point const &lpos) {

	bool const ret(lpos != last_lpos);
	last_lpos = lpos;
	return ret;
}

bool ground_mode_smap_data_t::needs_update(point const &lpos) {

	bool const has_dynamic(!tid || scene_smap_vbo_invalid || no_sparse_smap_update()); // Note: force two frames of updates the first time the smap is created by setting has_dynamic
	bool const ret(smap_data_t::needs_update(lpos) || has_dynamic || last_has_dynamic || voxel_shadows_updated); // Note: see view clipping in indexed_vntc_vect_t<T>::render()
	last_has_dynamic = has_dynamic;
	return ret;
}


void smap_data_t::create_shadow_map_for_light(int light, point const &lpos, cube_t const &bounds) { // Note: light is unused

	// setup render state
	bool const do_update(needs_update(lpos)); // must be called first, because this may indirectly update bounds
	xform_matrix const camera_mv_matrix(fgGetMVM()); // cache the camera modelview matrix before we change it
	fgPushMatrix();
	fgMatrixMode(FG_PROJECTION);
	fgPushMatrix();
	fgMatrixMode(FG_MODELVIEW);
	pdu = get_pt_cube_frustum_pdu(lpos, bounds, 1);
	texture_matrix = get_texture_matrix(camera_mv_matrix);
	check_gl_error(201);

	if (do_update) {
		// setup textures and framebuffer
		if (!tid) {
			bool const nearest(0); // nearest filter: sharper shadow edges, but needs more biasing
			setup_texture(tid, 0, 0, 0, 0, 0, nearest);
			set_shadow_tex_params();
			glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadow_map_sz, shadow_map_sz, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
		}
		// render from the light POV to a FBO, store depth values only
		enable_fbo(fbo_id, tid, 1);
		glViewport(0, 0, shadow_map_sz, shadow_map_sz);
		glClear(GL_DEPTH_BUFFER_BIT);
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
		// save state and update variables for fast rendering with correct clipping
		unsigned const orig_enabled_lights(enabled_lights);
		pos_dir_up const camera_pdu_(camera_pdu);
		camera_pdu     = pdu;
		enabled_lights = 0; // disable lighting so that shaders that auto-detect enabled lights don't try to do lighting
		ensure_filled_polygons();
		// render the scene
		check_gl_error(202);
		render_scene_shadow_pass(lpos);
		// restore state variables
		camera_pdu     = camera_pdu_;
		enabled_lights = orig_enabled_lights;
		disable_fbo();
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		set_standard_viewport();
		set_fill_mode();
	}
	restore_prev_mvm_pjm_state();
	// Now rendering from the camera POV, using the FBO to generate shadows
	check_gl_error(203);
}


void ground_mode_smap_data_t::render_scene_shadow_pass(point const &lpos) {

	point const camera_pos_(camera_pos);
	camera_pos = lpos;

	// add static objects
	if (coll_objects.drawn_ids.empty()) {
		// do nothing
	}
	else if (!smap_vertex_cache.vbo_valid()) {
		// only valid if drawing trees, small trees, and scenery separately
		vector<vert_wrap_t> verts;
		vector<pair<float, unsigned> > z_sorted;

		for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
			coll_obj const &c(coll_objects[*i]);
			assert(c.cp.draw);
			if (c.no_shadow_map()) continue;
			int ndiv(1);

			if (c.type == COLL_CUBE) {
				z_sorted.push_back(make_pair(-c.d[2][1], *i));
				continue;
			}
			else if (c.type == COLL_SPHERE) {
				ndiv = get_smap_ndiv(c.radius);
			}
			else if (c.type == COLL_CYLINDER || c.type == COLL_CYLINDER_ROT) {
				ndiv = get_smap_ndiv(max(c.radius, c.radius2));
			}
			c.get_shadow_triangle_verts(verts, ndiv);
		}
		smap_vertex_cache.end_block1(verts.size());
		sort(z_sorted.begin(), z_sorted.end());

		for (vector<pair<float, unsigned> >::const_iterator i = z_sorted.begin(); i != z_sorted.end(); ++i) {
			coll_objects[i->second].get_shadow_triangle_verts(verts, 1);
		}
		smap_vertex_cache.upload(verts);
	}
	smap_vertex_cache.render();
	render_models(1);
	render_voxel_data(1);

	if (!shadow_objs.empty()) { // add dynamic objects
		vector<vert_wrap_t> &dverts(smap_vertex_cache.dverts);
		shader_t shader;
		shader.set_vert_shader("vertex_xlate_scale");
		shader.set_frag_shader("color_only");
		shader.begin_shader();
		int const shader_loc(shader.get_uniform_loc("xlate_scale"));
		assert(shader_loc >= 0);
		bind_draw_sphere_vbo(0, 0); // no tex coords or normals

		for (vector<shadow_sphere>::const_iterator i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
			if (!pdu.sphere_visible_test(i->pos, i->radius)) continue;
			int const ndiv(get_smap_ndiv(i->radius));

			if (i->ctype != COLL_SPHERE) {
				assert((unsigned)i->cid < coll_objects.size());
				coll_objects[i->cid].get_shadow_triangle_verts(dverts, ndiv);
			}
			else {
				shader_t::set_uniform_vector4d(shader_loc, vector4d(i->pos, i->radius));
				draw_sphere_vbo_pre_bound(ndiv, 0);
			}
		}
		bind_vbo(0); // clear any bound sphere VBOs
		shader.set_uniform_vector4d(shader_loc, vector4d(all_zeros, 1.0)); // reset to identity transform
		smap_vertex_cache.render_dynamic();
		shader.end_shader();
	}
	// add snow, trees, scenery, and mesh
	if (snow_shadows) {draw_snow(1);} // slow
	draw_trees(1);
	draw_scenery(1, 1, 1);

	if ((display_mode & 0x01) && ground_effects_level != 0) { // draw mesh
		fgPushMatrix();
		float const val(1.0/dot_product(lpos.get_norm(), plus_z));
		fgTranslate(0.0, 0.0, -val*approx_pixel_width()); // translate down slightly to reduce shadow aliasing problems
		display_mesh(1);
		fgPopMatrix();
	}
	voxel_shadows_updated = 0;
	camera_pos = camera_pos_;
}


void create_shadow_map() {

	if (!shadow_map_enabled()) return; // disabled
	//RESET_TIME;

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
	orig_camera_pdu = camera_pdu;

	// set to shadow map state
	do_zoom  = 0;
	animate2 = 0; // disable any animations or generated effects
	display_mode &= ~0x08; // disable occlusion culling

	// check VBO
	if (scene_smap_vbo_invalid) {free_smap_vbo();}

	// render shadow maps to textures
	add_coll_shadow_objs();
	
	for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		point lpos;
		if (light_valid_and_enabled(l, lpos)) {smap_data[l].create_shadow_map_for_light(l, lpos, get_scene_bounds());}
	}
	scene_smap_vbo_invalid = 0;

	// restore old state
	check_gl_error(200);
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	//PRINT_TIME("Shadow Map Creation");
}


void free_shadow_map_textures() {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		smap_data[l].free_gl_state();
	}
	free_smap_vbo();
}




