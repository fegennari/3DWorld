// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "GL/glew.h"
#include "3DWorld.h"
#include "collision_detect.h" // for shadow_sphere
#include "gl_ext_arb.h"

using namespace std;

unsigned const SHADOW_MAP_SZ = 1024; // width/height - might need to be larger

int enable_shadow_maps(2); // 1 = dynamic shadows, 2 = dynamic + static shadows

extern bool have_drawn_cobj;
extern int window_width, window_height, animate2, display_mode;
extern vector<shadow_sphere> shadow_objs;
extern vector<coll_obj> coll_objects;


struct smap_data_t {
	unsigned tid, tu_id, fbo_id;
	float approx_pixel_width;
	pos_dir_up pdu;

	smap_data_t() : tid(0), tu_id(0), fbo_id(0), approx_pixel_width(0.0) {}

	int get_ndiv(float radius) {
		// FIXME: dynamic based on distance(camera, line(lpos, scene_center))?
		return min(N_SPHERE_DIV, max(3, int(radius/approx_pixel_width)));
	}
};

smap_data_t smap_data[NUM_LIGHT_SRC];


// ************ RENDER TO TEXTURE METHOD ************


unsigned get_shadow_map_tu_id(int light) {

	return smap_data[light].tu_id;
}

unsigned get_shadow_map_tid(int light) {

	return smap_data[light].tid;
}


void set_texture_matrix() {

	double modelView[16], projection[16];
	
	// This matrix transforms every coordinate x,y,z
	// x = x* 0.5 + 0.5 
	// y = y* 0.5 + 0.5 
	// z = z* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	const double bias[16] = {	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
		0.5, 0.5, 0.5, 1.0};
	
	// Grab modelview and projection matrices
	glGetDoublev(GL_MODELVIEW_MATRIX,  modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glLoadMatrixd(bias);
	
	// Concatating all matrice into one
	glMultMatrixd(projection);
	glMultMatrixd(modelView);
	
	// Go back to normal matrix mode
	glMatrixMode(GL_MODELVIEW);
}


void set_smap_shader_for_light(unsigned p, int light) {

	assert(light >= 0 && light < NUM_LIGHT_SRC);
	smap_data_t const &data(smap_data[light]);
	assert(data.tid > 0);
	add_uniform_int(p, append_array_ix(string("sm_tu_id"), light), data.tu_id);
	add_uniform_int(p, append_array_ix(string("sm_tex"),   light), data.tu_id);
	set_multitex(data.tu_id);
	bind_2d_texture(data.tid);
	set_multitex(0);
}


void set_smap_shader_for_all_lights(unsigned p) {

	point lpos; // unused

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
		if (!light_valid(0xFF, l, lpos)) continue;
		set_smap_shader_for_light(p, l);
	}
}


pos_dir_up get_light_pdu(point const &lpos, bool set_pers, bool do_look_at, float *sradius=NULL) {

	float const scene_z1(min(zbottom, czmin)), scene_z2(max(ztop, czmax)), scene_dz(scene_z2 - scene_z1);
	point const scene_center(0.0, 0.0, 0.5*(scene_z1 + scene_z2));
	float const scene_radius(sqrt(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE + scene_dz*scene_dz));
	cube_t const scene_bounds(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, scene_z1, scene_z2);
	point corners[8];
	get_cube_corners(scene_bounds.d, corners);
	float scene_radius2(0.0);

	for (unsigned i = 0; i < 8; ++i) {
		scene_radius2 = max(scene_radius2, pt_line_dist(corners[i], lpos, scene_center));
	}
	assert(scene_radius2 <= scene_radius);
	if (sradius) *sradius = scene_radius2;
	vector3d const light_dir((scene_center - lpos).get_norm()); // almost equal to lpos (point light)
	float const dist(p2p_dist(lpos, scene_center));
	vector3d up_dir(zero_vector);
	up_dir[get_min_dim(light_dir)] = 1.0;
	float const angle(atan2(scene_radius2, dist));
	pos_dir_up const pdu(lpos, light_dir, up_dir, tanf(angle)*SQRT2, sinf(angle), dist-scene_radius, dist+scene_radius, 1.0);

	if (set_pers) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(2.0*angle/TO_RADIANS, 1.0, pdu.near_, pdu.far_);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	if (do_look_at) {
		gluLookAt(lpos.x, lpos.y, lpos.z, scene_center.x, scene_center.y, scene_center.z, up_dir.x, up_dir.y, up_dir.z);
	}
	return pdu;
}


void draw_scene_bounds_and_light_frustum(point const &lpos) {

	glDisable(GL_LIGHTING);
	enable_blend();
	plus_z.do_glNormal(); // probably not needed

	// draw scene bounds
	glColor4f(1.0, 1.0, 1.0, 0.25);
	float const scene_z1(min(zbottom, czmin)), scene_z2(max(ztop, czmax));
	draw_simple_cube(cube_t(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, scene_z1, scene_z2), 0);

	// draw light frustum
	glColor4f(1.0, 1.0, 0.0, 0.25);
	get_light_pdu(lpos, 0, 0).draw_frustum();
	disable_blend();
	glEnable(GL_LIGHTING);
}


void create_shadow_fbo(unsigned &fbo_id, unsigned depth_tid) {
	
	// Create a framebuffer object
	glGenFramebuffers(1, &fbo_id);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);
	
	// Instruct openGL that we won't bind a color texture with the currently binded FBO
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	
	// Attach the texture to FBO depth attachment point
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_tid, 0);
	
	// Check FBO status
	GLenum const status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
	assert(status == GL_FRAMEBUFFER_COMPLETE);
	
	// Switch back to window-system-provided framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void create_shadow_map_for_light(int light, point const &lpos) {

	smap_data_t &data(smap_data[light]);
	data.tu_id = (6 + light); // Note: only 8 TUs guaranteed so we can have 2 lights

	// setup textures and framebuffer
	if (!data.tid) {
		setup_texture(data.tid, GL_MODULATE, 0, 0, 0, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_SZ, SHADOW_MAP_SZ, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
		glDisable(GL_TEXTURE_2D);
	}
	if (!data.fbo_id) create_shadow_fbo(data.fbo_id, data.tid);
	assert(data.fbo_id > 0);
	// Render from the light POV to a FBO, store depth values only
	glBindFramebuffer(GL_FRAMEBUFFER, data.fbo_id); // Rendering offscreen

	// setup render state
	glViewport(0, 0, SHADOW_MAP_SZ, SHADOW_MAP_SZ);
	glClear(GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	float sradius(0.0);
	camera_pos = lpos;
	data.pdu   = camera_pdu = get_light_pdu(lpos, 1, 1, &sradius);
	camera_pdu.valid = 0; // FIXME: should anything ever be out of the light view frustum? the camera when not over the mesh?
	data.approx_pixel_width = sradius / SHADOW_MAP_SZ;
	check_gl_error(201);

	// setup texture matrix
	set_multitex(data.tu_id);
	set_texture_matrix();
	disable_multitex_a();

	// render shadow geometry
	glDisable(GL_LIGHTING);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
	//glEnable(GL_CULL_FACE); glCullFace(GL_FRONT);
	WHITE.do_glColor();
	check_gl_error(202);

	if (enable_shadow_maps) { // add dynamic objects
		for (vector<shadow_sphere>::const_iterator i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
			if (!data.pdu.sphere_visible_test(i->pos, i->radius)) continue;
			int const ndiv(data.get_ndiv(i->radius));

			if (i->ctype != COLL_SPHERE) {
				assert((unsigned)i->cid < coll_objects.size());
				coll_objects[i->cid].simple_draw(ndiv);
			}
			else {
				// FIXME: use circle texture billboards
				draw_sphere_dlist(i->pos, i->radius, ndiv, 0);
			}
		}
	}
	if (enable_shadow_maps == 2) { // add static objects
		if (have_drawn_cobj) {
			int in_cur_prim(PRIM_UNSET);
			point center;
			float radius;

			for (vector<coll_obj>::const_iterator i = coll_objects.begin(); i != coll_objects.end(); ++i) {
				if (i->status != COLL_STATIC) continue;
				//if (i->no_draw()) continue;
				if (i->cp.color.alpha < MIN_SHADOW_ALPHA) continue;
				i->bounding_sphere(center, radius);
				if (!data.pdu.sphere_visible_test(center, radius)) continue;
				int ndiv(1);

				if (i->type == COLL_SPHERE) {
					ndiv = data.get_ndiv(radius);
				}
				else if (i->type == COLL_CYLINDER || i->type == COLL_CYLINDER_ROT) {
					ndiv = data.get_ndiv(max(i->radius, i->radius2));
				}
				in_cur_prim = i->simple_draw(ndiv, in_cur_prim, 1);
			}
			if (in_cur_prim >= 0) glEnd();
		}
		// FIXME: WRITE: render other static objects such as mesh
		// FIXME: remember to handle trasparency
	}
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	//glDisable(GL_CULL_FACE); glCullFace(GL_BACK);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glEnable(GL_LIGHTING);

	// Now rendering from the camera POV, using the FBO to generate shadows
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	check_gl_error(203);
}


void create_shadow_map() {

	if (!enable_shadow_maps) return; // disabled
	//RESET_TIME;

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
	point const camera_pos_(camera_pos);
	pos_dir_up const camera_pdu_(camera_pdu);

	// set to shadow map state
	do_zoom       = 0;
	animate2      = 0; // disable any animations or generated effects
	display_mode &= ~0x08; // disable occlusion culling

	// render shadow maps to textures
	point lpos;
	bool smap_used(0);
	
	for (unsigned is_dynamic = 0; is_dynamic < 2; ++is_dynamic) { // {static, dynamic}
		for (int l = 0; l < NUM_LIGHT_SRC; ++l) { // {sun, moon}
			if (!light_valid(0xFF, l, lpos)) continue;
			create_shadow_map_for_light(l, lpos);
			smap_used = 1;
		}
	}

	// restore old state
	if (smap_used) {
		glDisable(GL_TEXTURE_2D);
		glViewport(0, 0, window_width, window_height);
		glClear(GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	}
	check_gl_error(200);
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	camera_pos   = camera_pos_;
	camera_pdu   = camera_pdu_;
	//PRINT_TIME("Shadow Map Creation");
}


void free_shadow_map_textures() {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		free_texture(smap_data[l].tid);
		
		if (smap_data[l].fbo_id > 0) {
			glDeleteFramebuffers(1, &smap_data[l].fbo_id);
			smap_data[l].fbo_id = 0;
		}
	}
}




