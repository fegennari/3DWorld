// 3D World - GL EXT/ARB extension interface code
// by Frank Gennari
// 4/25/02
#include "3DWorld.h"
#include "gl_ext_arb.h"


void init_glew() {

// MacOSX check here, placeholder for eventual cross-platform porting
#if ((defined(__MACH__))&&(defined(__APPLE__)))
	// nothing
#else
	GLenum const err(glewInit());

	if (GLEW_OK != err) {
	  fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	  assert(0);
	}
#endif
	glEnable(GL_MULTISAMPLE); // only works when using a multisampling graphics context
}


// ***************** MULTITEXTURING *****************

unsigned const MAX_MULTITEX = 32; // max is GL_TEXTURE31

unsigned max_used_multitex(0);
bool multitex_enabled[MAX_MULTITEX] = {0};


void set_multitex(unsigned tu_id) {

	static bool inited(0);
	if (!inited) setup_multitexture();
	inited = 1;
	assert(tu_id < MAX_MULTITEX); // Note: Assumes textures are defined sequentially
	multitex_enabled[tu_id] = 1;
	max_used_multitex = max(max_used_multitex, (tu_id+1));
	glActiveTexture(GL_TEXTURE0 + tu_id);
}


void select_multitex(int id, unsigned tu_id, bool enable) {

	set_multitex(tu_id);
	select_texture(id, enable);
}


void disable_multitex(unsigned tu_id, bool reset) {

	assert(tu_id < max_used_multitex);
	multitex_enabled[tu_id] = 0;
	if ((tu_id+1) == max_used_multitex) --max_used_multitex;
	glActiveTexture((GL_TEXTURE0 + tu_id));
	glDisable(GL_TEXTURE_2D);
	if (reset) glActiveTexture(GL_TEXTURE0); // end back at texture 0
}


void disable_multitex_a() {

	for (unsigned i = 0; i < max_used_multitex; ++i) {
		if (multitex_enabled[i]) disable_multitex(i, 0);
	}
	glActiveTexture(GL_TEXTURE0); // end back at texture 0
	max_used_multitex = 0;
}


void multitex_coord_n(unsigned tu_id, float const *v, unsigned num) {

	unsigned const id(GL_TEXTURE0 + tu_id);
	switch (num) {
		case 1: glMultiTexCoord1fv(id, v); break;
		case 2: glMultiTexCoord2fv(id, v); break;
		case 3: glMultiTexCoord3fv(id, v); break;
		case 4: glMultiTexCoord4fv(id, v); break;
		default: assert(0);
	}
}


void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id) {

	assert(tu_id < max_used_multitex);
	glMultiTexCoord2f((GL_TEXTURE0 + tu_id), s, t);
}


void multitex_coord2f_a(GLfloat s, GLfloat t) {

	for (unsigned i = 0; i < max_used_multitex; ++i) {
		if (multitex_enabled[i]) multitex_coord2f(s, t, i);
	}
}


void setup_multitexture() { // Windows specific

	if (!glMultiTexCoord2f) {
		cout << "*** Can't find GL_multitexture extension ***" << endl;
		assert(0);
	}
	//glClientActiveTexture(GL_TEXTURE0);
	//glTexCoordPointer(2, GL_FLOAT, 0, tp0);
}


// ***************** 3D TEXTURES *****************


void bind_3d_texture(unsigned tid) {

	glBindTexture(GL_TEXTURE_3D, tid);
	assert(glIsTexture(tid));
}


GLenum get_format(unsigned ncomp) {
	return ((ncomp == 1) ? GL_LUMINANCE : ((ncomp == 4) ? GL_RGBA : GL_RGB));
}


unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data, int filter, int wrap) {

	assert(data.size() == ncomp*xsz*ysz*zsz);
	unsigned tid(0);
	glGenTextures(1, &tid);
	bind_3d_texture(tid);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, filter);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, filter);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, wrap);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, wrap);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, wrap);
	glTexImage3D(GL_TEXTURE_3D, 0, ncomp, xsz, ysz, zsz, 0, get_format(ncomp), GL_UNSIGNED_BYTE, &data.front());
	return tid;
}


void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data)
{
	assert(glIsTexture(tid));
	bind_3d_texture(tid);
	glTexSubImage3D(GL_TEXTURE_3D, 0, xoff, yoff, zoff, xsz, ysz, zsz, get_format(ncomp), GL_UNSIGNED_BYTE, data);
}


// ***************** FOG_COORD *****************


void setup_fog_coord() {

	if (!glFogCoordf) {
		cout << "*** Can't find fog_coord extension ***" << endl;
		assert(0);
	}
}


void set_fog_coord(GLfloat val) {
	glFogCoordf(val);
}

void enable_fog_coord() {
	glFogi(GL_FOG_COORDINATE_SOURCE, GL_FOG_COORDINATE);
}

void disable_fog_coord() {
	glFogi(GL_FOG_COORDINATE_SOURCE, GL_FRAGMENT_DEPTH); // is this correct?
}


// ***************** VBOs *****************


bool setup_gen_buffers() {

	static int retval(2);
	
	if (retval == 2) { // not yet setup
		if (!glGenBuffers) {
			cout << "*** glGenBuffers is not supported ***" << endl;
			retval = 0;
		}
		else {
			retval = 1;
		}
	}
	return (retval != 0);
}

unsigned create_vbo() {
	unsigned id;
	glGenBuffers(1, &id);
	assert(id > 0);
	return id;
}

int get_buffer_target(bool is_index) {return (is_index ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER);}

void bind_vbo(unsigned id, bool is_index) {
	glBindBuffer(get_buffer_target(is_index), id);
}

void delete_vbo(unsigned id) {
	if (id == 0) return;
	glDeleteBuffers(1, &id);
}

void upload_vbo_data(void const *const data, size_t size, bool is_index) {
	// hard coded for drawing of static data
	glBufferData(get_buffer_target(is_index), size, data, GL_STATIC_DRAW);
}

void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index) {
	glBufferSubData(get_buffer_target(is_index), offset, size, data);
}


// ***************** FBOs *****************


void create_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo) {
	
	// Create a framebuffer object
	glGenFramebuffers(1, &fbo_id);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);
	
	if (is_depth_fbo) {
		// Instruct openGL that we won't bind a color texture with the currently binded FBO
		glDrawBuffer(GL_NONE);
		glReadBuffer(GL_NONE);
	}
	
	// Attach the texture to FBO depth or color attachment point
	glFramebufferTexture2D(GL_FRAMEBUFFER, (is_depth_fbo ? GL_DEPTH_ATTACHMENT : GL_COLOR_ATTACHMENT0), GL_TEXTURE_2D, tid, 0);
	
	// Check FBO status
	GLenum const status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
	assert(status == GL_FRAMEBUFFER_COMPLETE);
	
	// Switch back to window-system-provided framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void enable_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo) {

	assert(glIsTexture(tid));
	if (!fbo_id) create_fbo(fbo_id, tid, is_depth_fbo);
	assert(fbo_id > 0);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_id); // Rendering offscreen
}


void disable_fbo() {

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void free_fbo(unsigned &fbo_id) {

	if (fbo_id > 0) glDeleteFramebuffers(1, &fbo_id);
	fbo_id = 0;
}


// ***************** Other *****************


bool gen_mipmaps() {

	if (!glGenerateMipmap) return 0;
	glGenerateMipmap(GL_TEXTURE_2D);
	return 1;
}


