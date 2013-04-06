// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "3DWorld.h"
#include "collision_detect.h"
#include "gl_ext_arb.h"
#include "transform_obj.h" // for xform_matrix
#include "shaders.h"
#include "model3d.h"


bool const ENABLE_DLIST = 1;

bool scene_dlist_invalid(0);
unsigned shadow_map_sz(0), smap_dlist(0);
pos_dir_up orig_camera_pdu;

extern int window_width, window_height, animate2, display_mode, ground_effects_level;
extern vector<shadow_sphere> shadow_objs;
extern coll_obj_group coll_objects;

void draw_trees(bool shadow_only=0);


struct smap_data_t {
	unsigned tid, tu_id, fbo_id;
	pos_dir_up pdu;

	smap_data_t() : tid(0), tu_id(0), fbo_id(0) {}

	void free_gl_state() {
		free_texture(tid);
		free_fbo(fbo_id);
	}
	void create_shadow_map_for_light(int light, point const &lpos);
};

smap_data_t smap_data[NUM_LIGHT_SRC];


bool shadow_map_enabled() {
	return (shadow_map_sz > 0);
}

unsigned get_shadow_map_tu_id(int light) {
	return smap_data[light].tu_id;
}

unsigned get_shadow_map_tid(int light) {
	return smap_data[light].tid;
}

float approx_pixel_width() {
	return 0.5*sqrt(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE) / shadow_map_sz;
}

int get_smap_ndiv(float radius) {
	// dynamic based on distance(camera, line(lpos, scene_center))?
	return min(N_SPHERE_DIV, max(3, int(radius/approx_pixel_width())));
}


void free_smap_dlist() {

	if (glIsList(smap_dlist)) glDeleteLists(smap_dlist, 1);
	smap_dlist = 0;
}


struct matrix4x4d : public xform_matrix {

	struct matrix3x3d {
		vector3d_d x,y,z;

		double get_determinant() const {
			return (x.x*(y.y*z.z - y.z*z.y) - y.x*(x.y*z.z - x.z*z.y) + z.x*(x.y*y.z - x.z*y.y));
		}
	};

	matrix3x3d get_sub_matrix(unsigned x, unsigned y) const {
		matrix3x3d tmp;
		unsigned xoffset(0);

		for (unsigned i = 0; i < 4; i++) {
			if (i == x) continue;
			unsigned yoffset(0);

			for (unsigned j = 0; j < 4; j++) {
				if (j == y) continue;
				*(((double*) &tmp) + xoffset*3 + yoffset) = m[i*4 + j];
				yoffset++;
			}
			xoffset++;
		}
		return tmp;
	}

	double get_determinant() const {
		double result(0.0), i(1.0);

		for (unsigned n = 0; n < 4; n++, i *= -1.0) {
			result += m[n] * get_sub_matrix(0, n).get_determinant() * i;
		}
		return result;
	}

	matrix4x4d get_inverse() {
		matrix4x4d inverse;
		double const m4determinant(get_determinant());
		assert(fabs(m4determinant) > TOLERANCE);

		for (unsigned i = 0; i < 4; i++) {
			for (unsigned j = 0; j < 4; j++) {
				int const sign(1-((i+j)&1)*2);
				inverse.m[i + j*4] = get_sub_matrix(i, j).get_determinant() * sign / m4determinant;
			}
		}
		return inverse;
	}
};


void set_texture_matrix(matrix4x4d &camera_mv_matrix) {

	matrix4x4d modelView, projection;
	
	// This matrix transforms every coordinate {x,y,z} to {x,y,z}* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	const double bias[16] = {	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
		0.5, 0.5, 0.5, 1.0};
	
	// Grab modelview and projection matrices
	modelView.assign_mv_from_gl();
	projection.assign_pj_from_gl();
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glLoadMatrixd(bias);
	
	// Concatating all matrice into one
	projection.apply();
	modelView.apply();
	camera_mv_matrix.get_inverse().apply();
	
	// Go back to normal matrix mode
	glMatrixMode(GL_MODELVIEW);
}


bool set_smap_shader_for_light(shader_t &s, int light, float z_bias) {

	if (!shadow_map_enabled()) return 0;
	assert(light >= 0 && light < NUM_LIGHT_SRC);
	if (!glIsEnabled(GL_LIGHT0 + light)) return 0;
	point lpos; // unused
	smap_data_t const &data(smap_data[light]);
	assert(data.tid > 0);
	s.add_uniform_float("z_bias", z_bias);
	s.add_uniform_int  (append_ix(string("sm_tu_id"), light, 0), data.tu_id);
	s.add_uniform_int  (append_ix(string("sm_tex"),   light, 0), data.tu_id);
	s.add_uniform_float(append_ix(string("sm_scale"), light, 0), (light_valid(0xFF, light, lpos) ? 1.0 : 0.0));
	set_active_texture(data.tu_id);
	bind_2d_texture(data.tid);
	set_active_texture(0);
	return 1;
}


void set_smap_shader_for_all_lights(shader_t &s, float z_bias) {

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		set_smap_shader_for_light(s, l, z_bias);
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

#if 1 // tighter bounds / higher quality / slower
	vector3d dirs[2]; // x, y (up)
	orthogonalize_dir(up_dir, light_dir, dirs[1], 1);
	cross_product(light_dir, dirs[1], dirs[0]);
	float rx(0.0), ry(0.0);

	for (unsigned i = 0; i < 8; ++i) {
		rx = max(rx, fabs(dot_product(dirs[0], (corners[i] - pos))));
		ry = max(ry, fabs(dot_product(dirs[1], (corners[i] - pos))));
	}
	float const angle(atan2(ry, dist)), aspect(rx/ry);
#else
	float radius2(0.0);

	for (unsigned i = 0; i < 8; ++i) {
		radius2 = max(radius2, pt_line_dist(corners[i], pos, center));
	}
	assert(radius2 <= radius); // can fail for cloud frustum?
	float const angle(atan2(radius2, dist)), aspect(1.0);
#endif
	pos_dir_up const pdu(pos, light_dir, up_dir, tanf(angle)*SQRT2, sinf(angle), max(NEAR_CLIP, dist-radius), dist+radius, aspect);

	if (set_matrix) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(2.0*angle/TO_RADIANS, aspect, pdu.near_, pdu.far_);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(pos.x, pos.y, pos.z, center.x, center.y, center.z, up_dir.x, up_dir.y, up_dir.z);
	}
	return pdu;
}


cube_t get_scene_bounds() {

	return cube_t(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, min(zbottom, czmin), max(ztop, czmax));
}


void draw_scene_bounds_and_light_frustum(point const &lpos) {

	glDisable(GL_LIGHTING);
	enable_blend();

	// draw scene bounds
	glColor4f(1.0, 1.0, 1.0, 0.25);
	draw_simple_cube(get_scene_bounds(), 0);

	// draw light frustum
	glColor4f(1.0, 1.0, 0.0, 0.25);
	get_pt_cube_frustum_pdu(lpos, get_scene_bounds(), 0).draw_frustum();
	disable_blend();
	glEnable(GL_LIGHTING);
}


void set_shadow_tex_params() {

	// This is to allow usage of shadow2DProj function in the shader
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
}


void smap_data_t::create_shadow_map_for_light(int light, point const &lpos) {

	tu_id = (6 + light); // Note: currently used with 2 lights, up to TU7

	// setup textures and framebuffer
	if (!tid) {
		bool const nearest(0); // nearest filter: sharper shadow edges, but needs more biasing
		setup_texture(tid, GL_MODULATE, 0, 0, 0, 0, 0, nearest);
		set_shadow_tex_params();
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadow_map_sz, shadow_map_sz, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
	}

	// Render from the light POV to a FBO, store depth values only
	enable_fbo(fbo_id, tid, 1);

	// setup render state
	matrix4x4d camera_mv_matrix;
	camera_mv_matrix.assign_mv_from_gl(); // cache the camera modelview matrix before we change it
	glViewport(0, 0, shadow_map_sz, shadow_map_sz);
	glClear(GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	float sradius(0.0);
	camera_pos = lpos;
	pdu        = camera_pdu = get_pt_cube_frustum_pdu(lpos, get_scene_bounds(), 1);
	check_gl_error(201);

	// setup texture matrix
	set_active_texture(tu_id);
	set_texture_matrix(camera_mv_matrix);
	set_active_texture(0);

	// render shadow geometry
	glDisable(GL_LIGHTING);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
	WHITE.do_glColor();
	check_gl_error(202);

	// add static objects
	if (coll_objects.drawn_ids.empty()) {
		// do nothing
	}
	else if (smap_dlist) {
		assert(glIsList(smap_dlist));
		glCallList(smap_dlist);
	}
	else {
		if (ENABLE_DLIST) {
			smap_dlist = glGenLists(1);
			glNewList(smap_dlist, GL_COMPILE_AND_EXECUTE);
		}
		glEnable(GL_CULL_FACE);

		for (unsigned n = 0; n < 2; ++n) { // {cubes+culling, others+no cuqlling}
			int in_cur_prim(PRIM_UNSET);

			// FIXME: sort cubes in decreasing z-value?
			// only valid if drawing trees, small trees, and scenery separately
			for (cobj_id_set_t::const_iterator i = coll_objects.drawn_ids.begin(); i != coll_objects.drawn_ids.end(); ++i) {
				coll_obj const &c(coll_objects[*i]);
				assert(c.cp.draw);
				if ((c.type == COLL_CUBE) == n)   continue;
				if (c.no_draw() || c.no_shadow()) continue;
				int ndiv(1);

				if (c.type == COLL_SPHERE) {
					ndiv = get_smap_ndiv(c.radius);
				}
				else if (c.type == COLL_CYLINDER || c.type == COLL_CYLINDER_ROT) {
					ndiv = get_smap_ndiv(max(c.radius, c.radius2));
				}
				in_cur_prim = c.simple_draw(ndiv, in_cur_prim, 1, ENABLE_DLIST);
			}
			if (in_cur_prim >= 0) glEnd();
			glDisable(GL_CULL_FACE);
		} // for n
		if (ENABLE_DLIST) glEndList();
	}
	render_models(1);
	render_voxel_data(1);

	// add dynamic objects
	for (vector<shadow_sphere>::const_iterator i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
		if (!pdu.sphere_visible_test(i->pos, i->radius)) continue;
		int const ndiv(get_smap_ndiv(i->radius));

		if (i->ctype != COLL_SPHERE) {
			assert((unsigned)i->cid < coll_objects.size());
			coll_objects[i->cid].simple_draw(ndiv, PRIM_DISABLED, 1, 0);
		}
		else {
			draw_sphere_dlist(i->pos, i->radius, ndiv, 0); // use circle texture billboards?
		}
	}
	draw_trees(1);
	draw_scenery(1, 1, 1);

	if ((display_mode & 0x01) && ground_effects_level != 0) { // draw mesh
		glPushMatrix();
		float const val(1.0/dot_product(lpos.get_norm(), plus_z));
		glTranslatef(0.0, 0.0, -val*approx_pixel_width()); // translate down slightly to reduce shadow aliasing problems
		display_mesh(1);
		glPopMatrix();
	}
	
	// reset state
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glEnable(GL_LIGHTING);

	// Now rendering from the camera POV, using the FBO to generate shadows
	disable_fbo();
	check_gl_error(203);
}


void create_shadow_map() {

	if (!shadow_map_enabled()) return; // disabled
	//RESET_TIME;

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
	point const camera_pos_(camera_pos);
	orig_camera_pdu = camera_pdu;

	// set to shadow map state
	do_zoom  = 0;
	animate2 = 0; // disable any animations or generated effects
	display_mode &= ~(0x08 | 0x0100); // disable occlusion culling and leaf wind

	// check dlist
	if (scene_dlist_invalid) {
		free_smap_dlist();
		scene_dlist_invalid = 0;
	}

	// render shadow maps to textures
	add_coll_shadow_objs();
	point lpos;
	
	for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		if (!light_valid(0xFF, l, lpos) || !glIsEnabled(GL_LIGHT0 + l)) continue;
		smap_data[l].create_shadow_map_for_light(l, lpos);
	}

	// restore old state
	set_standard_viewport();
	check_gl_error(200);
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	camera_pos   = camera_pos_;
	camera_pdu   = orig_camera_pdu;
	//PRINT_TIME("Shadow Map Creation");
}


void free_shadow_map_textures() {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		smap_data[l].free_gl_state();
	}
	free_smap_dlist();
}




