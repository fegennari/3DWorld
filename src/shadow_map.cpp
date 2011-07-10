// 3D World - Shadow Mapping using Shaders
// by Frank Gennari
// 1/21/11
#include "GL/glew.h"
#include "3DWorld.h"
#include "collision_detect.h" // for shadow_sphere
#include "gl_ext_arb.h"

using namespace std;

unsigned const SHADOW_MAP_SZ = 1024; // width/height - might need to be larger

unsigned fbo_id(0), depth_tid(0);

extern int window_width, window_height, animate2, display_mode;
extern vector<shadow_sphere> shadow_objs;


// ************ RENDER TO TEXTURE METHOD ************


unsigned create_shadow_map_for_light(int light, point const &lpos, unsigned tu_id, bool is_dynamic) { // return a tid

	assert(tu_id < 8); // assuming we have 8 texture units
	float const scene_z1(min(zbottom, czmin)), scene_z2(max(ztop, czmax));
	float const scene_radius(max(max(X_SCENE_SIZE, Y_SCENE_SIZE), 0.5f*(scene_z2 - scene_z1)));
	point const scene_center(0.0, 0.0, 0.5*(scene_z1 + scene_z2));
	vector3d const light_dir(scene_center - lpos); // almost equal to lpos (point light)
	float const dist(p2p_dist(lpos, scene_center));

	float const angle(0.5*TO_RADIANS*PERSP_ANGLE); // FIXME: what is the angle
	float const A((float)window_width/(float)window_height), diag(sqrt(1.0 + A*A)); // aspect ratio
	camera_pdu = pos_dir_up(lpos, light_dir, plus_z, tanf(angle)*diag, sinf(angle), dist-scene_radius, dist+scene_radius);
	camera_pdu.valid = 0; // FIXME: should anything ever be out of the light view frustum? the camera when not over the mesh?

	// FIXME: WRITE: setup render state

	if (is_dynamic) {
		//if (shadow_objs.empty()) {} // is this check necessary?
			
		for (vector<shadow_sphere>::const_iterator i = shadow_objs.begin(); i != shadow_objs.end(); ++i) {
			int const ndiv(N_SPHERE_DIV); // FIXME: dynamic based on distance(camera, line(lpos, scene_center))?
			draw_sphere_dlist(i->pos, i->radius, ndiv, 0);
		}
	}
	else {
		// FIXME: WRITE: render static objects (cobjs, etc)
		draw_coll_surfaces(1, 0); // solid, no transparent?
	}
	unsigned tid(0);
	// FIXME: WRITE: create textures
	return tid;
}


void create_shadow_map(bool is_dynamic) {

	// Note 1: We likely want both static and dynamic/per frame shadow maps
	// Note 2: We probably want to support at least two light sources, sun and moon
	// Note 3: There are 2-3 texture units available: 6, 7, and maybe 5 if the lflow texture is unused
	// We can't really support 4 shadow maps with 3 TUs, so something has to be done here

	// The plan:
	// 1. Start with dynamic shadows for grass with one light source
	// 2. Add cobj dynamic shadows
	// 3. Add other light sources (Sun + Moon)
	// 4. Add static light sources for a single light
	// 5. Resolve the issues above for static + dynamic + multiple lights

	// save state
	int const do_zoom_(do_zoom), animate2_(animate2), display_mode_(display_mode);
	pos_dir_up const camera_pdu_(camera_pdu);

	// set to shadow map state
	do_zoom       = 0;
	animate2      = 0; // disable any animations or generated effects
	display_mode &= ~0x08; // disable occlusion culling

	// render shadow maps to textures
	point lpos;

	for (int l = 0; l < NUM_LIGHT_SRC; ++l) {
		if (!light_valid(0xFF, l, lpos)) continue;
		unsigned const shadow_tid(create_shadow_map_for_light(l, lpos, 6+l, is_dynamic)); // FIXME: what about is_dynamic?
		// FIXME: use shadow_tid
	}

	// restore old state
	do_zoom      = do_zoom_;
	animate2     = animate2_;
	display_mode = display_mode_;
	camera_pdu   = camera_pdu_;
}


// ************ FBO METHOD ************


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


void set_texture_matrix() {

	// FIXME: GL transforms will invalidate the texture matrix
	static double modelView[16];
	static double projection[16];
	
	// This is matrix transform every coordinate x,y,z
	// x = x* 0.5 + 0.5 
	// y = y* 0.5 + 0.5 
	// z = z* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	const GLdouble bias[16] = {	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
		0.5, 0.5, 0.5, 1.0};
	
	// Grab modelview and transformation matrices
	glGetDoublev(GL_MODELVIEW_MATRIX, modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glMatrixMode(GL_TEXTURE);
	set_multitex(7);
	glLoadIdentity();
	glLoadMatrixd(bias);
	
	// Concatating all matrice into one
	glMultMatrixd(projection);
	glMultMatrixd(modelView);
	
	// Go back to normal matrix mode
	glMatrixMode(GL_MODELVIEW);
}


void setup_matrices(point const &pos, point const &look_at) {

	set_perspective(PERSP_ANGLE, 1.0);
	gluLookAt(pos.x, pos.y, pos.z, look_at.x, look_at.y, look_at.z, 0, 1, 0);
}


void draw_scene() {

	// placeholder for actual scene drawing
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
	draw_scene();
	
	// Save modelview/projection matrice into texture7, also add a biais
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



