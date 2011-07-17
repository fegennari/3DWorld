// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "GL/glew.h"
#include "3DWorld.h"
#include "collision_detect.h" // for shadow_sphere
#include "gl_ext_arb.h"

using namespace std;

unsigned const SHADOW_MAP_SZ = 1024; // width/height - might need to be larger

int enable_shadow_maps(1); // 1 = dynamic shadows, 2 = dynamic + static shadows
unsigned fbo_id(0), depth_tid(0);

extern bool have_drawn_cobj;
extern int window_width, window_height, animate2, display_mode;
extern vector<shadow_sphere> shadow_objs;
extern vector<coll_obj> coll_objects;


struct smap_data_t {
	bool last_no_dynamic;
	unsigned tid, tu_id, smap_sz;
	pos_dir_up pdu;

	smap_data_t() : last_no_dynamic(0), tid(0), tu_id(0), smap_sz(SHADOW_MAP_SZ) {}
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


pos_dir_up get_light_pdu(point const &lpos, bool set_pers, bool do_look_at) {

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


void create_shadow_map_for_light(int light, point const &lpos) {

	bool no_dynamic(shadow_objs.empty());
	smap_data_t &data(smap_data[light]);
	// FIXME: need to invalidate when static cobjs change
	if (data.tid > 0 && no_dynamic && data.last_no_dynamic) return; // no change to shadow map
	data.tu_id = (6 + light); // Note: only 8 TUs guaranteed so we can have 2 lights

	// determine shadow map size
	if (!data.tid) {
		unsigned const max_viewport(min(window_width, window_height));
		data.smap_sz = SHADOW_MAP_SZ;

		if (data.smap_sz > max_viewport) {
			cout << "Warning: Using smaller shadow map size of " << max_viewport << " due to small render window" << endl;
			data.smap_sz = max_viewport;
		}
	}

	// setup render state
	glViewport(0, 0, data.smap_sz, data.smap_sz);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	data.pdu = camera_pdu = get_light_pdu(lpos, 1, 1);
	camera_pdu.valid = 0; // FIXME: should anything ever be out of the light view frustum? the camera when not over the mesh?
	check_gl_error(201);

	// setup texture matrix
	set_multitex(data.tu_id);
	set_texture_matrix();
	disable_multitex_a();

	// render shadow geometry
	no_dynamic = 1;
	glDisable(GL_LIGHTING);
	WHITE.do_glColor();

	if (enable_shadow_maps) { // add dynamic objects
		for (vector<shadow_sphere>::const_iterator i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
			if (!data.pdu.sphere_visible_test(i->pos, i->radius)) continue;
			int const ndiv(N_SPHERE_DIV); // FIXME: dynamic based on distance(camera, line(lpos, scene_center))?

			if (i->ctype != COLL_SPHERE) {
				assert((unsigned)i->cid < coll_objects.size());
				coll_objects[i->cid].simple_draw(ndiv);
			}
			else {
				// FIXME: use circle texture billboards
				draw_sphere_dlist(i->pos, i->radius, ndiv, 0);
			}
			no_dynamic = 0;
		}
	}
	if (enable_shadow_maps == 2) { // add static objects
		if (have_drawn_cobj) {
			point center;
			float radius;

			for (unsigned i = 0; i < coll_objects.size(); ++i) {
				if (coll_objects[i].no_draw()) continue;
				if (coll_objects[i].cp.color.alpha < MIN_SHADOW_ALPHA) continue;
				coll_objects[i].bounding_sphere(center, radius);
				if (!data.pdu.sphere_visible_test(center, radius)) continue;
				int const ndiv(N_SPHERE_DIV); // is this enough/too many?
				coll_objects[i].simple_draw(ndiv);
			}
		}
		// FIXME: WRITE: render other static objects such as mesh
		// FIXME: remember to handle trasparency
	}
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_LIGHTING);
	check_gl_error(202);

	// setup textures
	if (!data.tid) {
		setup_texture(data.tid, GL_MODULATE, 0, 0, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, data.smap_sz, data.smap_sz, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
	}
	else {
		if (no_dynamic && data.last_no_dynamic) return; // no change
		bind_2d_texture(data.tid);
	}
	glReadBuffer(GL_BACK);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, data.smap_sz, data.smap_sz);
	check_gl_error(203);
	data.last_no_dynamic = no_dynamic;
}


void create_shadow_map() {

	if (!enable_shadow_maps) return; // disabled
	//RESET_TIME;

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
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
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	}
	check_gl_error(200);
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	camera_pdu   = camera_pdu_;
	//PRINT_TIME("Shadow Map Creation");
}


void free_shadow_map_textures() {

	for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
		free_texture(smap_data[l].tid);
	}
}


// ************ FBO METHOD - INCOMPLETE ************


// Note: Can be done similar to inf terrain mode water reflections
void create_shadow_fbo() {
	
	// Try to use a texture depth component
	assert(depth_tid == 0);
	glGenTextures(1, &depth_tid);
	glBindTexture(GL_TEXTURE_2D, depth_tid);
	
	// GL_LINEAR does not make sense for the depth texture.
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	
	// Remove artifact on the edges of the shadowmap
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	//glTexParameterfv( GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor );
	
	// No need to force GL_DEPTH_COMPONENT24, drivers usually give you the max precision if available 
	glTexImage2D( GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_SZ, SHADOW_MAP_SZ, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	
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


void setup_shadow_fbo() {

	// FIXME: need to reset fbo_id and depth_tid on minimize/maximize
	if (fbo_id == 0) create_shadow_fbo();
	assert(fbo_id    > 0);
	assert(depth_tid > 0);

	// This is important, if not here, the FBO's depthbuffer won't be populated.
	glEnable(GL_DEPTH_TEST);
	glClearColor(0,0,0,1.0f);
	glEnable(GL_CULL_FACE);
}


void setup_matrices(point const &pos, point const &look_at) {

	set_perspective(PERSP_ANGLE, 1.0);
	gluLookAt(pos.x, pos.y, pos.z, look_at.x, look_at.y, look_at.z, 0, 1, 0);
}


void render_to_shadow_fbo(point const &lpos) {

	setup_shadow_fbo();
	
	// First step: Render from the light POV to a FBO, store depth values only
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_id); // Rendering offscreen
	
	// Use the fixed pipeline to render to the depthbuffer
	// In the case we render the shadowmap to a higher resolution, the viewport must be modified accordingly.
	glViewport(0, 0, SHADOW_MAP_SZ, SHADOW_MAP_SZ);
	
	// Clear previous frame values
	glClear(GL_DEPTH_BUFFER_BIT);
	
	// Disable color rendering, we only want to write to the Z-Buffer
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); 
	
	setup_matrices(lpos, camera_origin);
	
	// Culling switching, rendering only backface, this is done to avoid self-shadowing
	glCullFace(GL_FRONT);
	//draw_scene();
	
	// Save modelview/projection matrice into texture7, also add a biais
	set_multitex(7);
	set_texture_matrix();
	
	// Now rendering from the camera POV, using the FBO to generate shadows
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	
	glViewport(0, 0, window_width, window_height);
	
	// Enable color write (previously disabled for light POV z-buffer rendering)
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE); 
	
	// Clear previous frame values
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Use the shadow shader
	//set_shader_prog("", ""); // FIXME
	set_multitex(7); // using depth texture 7
	glBindTexture(GL_TEXTURE_2D, depth_tid);
	glCullFace(GL_BACK);
	set_multitex(0);
	//unset_shader_prog();

	// Back to camera space
	
	// DEBUG only. this piece of code draws the depth buffer onscreen
#if 1
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-window_width/2, window_width/2, -window_height/2, window_height/2, 1, 20);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glColor4f(1,1,1,1);
	set_multitex(0);
	glBindTexture(GL_TEXTURE_2D, depth_tid);
	glEnable(GL_TEXTURE_2D);
	draw_tquad(window_width/2, window_height/2, -1, 1);
	glDisable(GL_TEXTURE_2D);
#endif

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	check_zoom();
	glDisable(GL_CULL_FACE);
	// draw the rest of the scene
}



