// 3D World - GL EXT/ARB extension interface code
// by Frank Gennari
// 4/25/02
#include "GL/glew.h" // must be included first
#include "3DWorld.h"
#include "gl_ext_arb.h"


void init_glew() {

	GLenum const err(glewInit());

	if (GLEW_OK != err) {
	  fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	  assert(0);
	}
	glEnable(GL_MULTISAMPLE);
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


void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id) {

	assert(tu_id < max_used_multitex);
	glMultiTexCoord2f((GL_TEXTURE0 + tu_id), s, t);
	//glTexCoord2f(s, t);
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
	//glEnableClientState(GL_TEXTURE_COORD_ARRAY);
}


// ***************** 3D TEXTURES *****************


unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data) {

	assert(data.size() == ncomp*xsz*ysz*zsz);
	unsigned tid(0);
	glGenTextures(1, &tid);
	glBindTexture(GL_TEXTURE_3D, tid);
	assert(glIsTexture(tid));
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	GLenum const format((ncomp == 4) ? GL_RGBA : GL_RGB);
	glTexImage3D(GL_TEXTURE_3D, 0, ncomp, xsz, ysz, zsz, 0, format, GL_UNSIGNED_BYTE, &data.front());
	return tid;
}


void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data)
{
	glBindTexture(GL_TEXTURE_3D, tid);
	assert(glIsTexture(tid));
	GLenum const format((ncomp == 4) ? GL_RGBA : GL_RGB);
	glTexSubImage3D(GL_TEXTURE_3D, 0, xoff, yoff, zoff, xsz, ysz, zsz, format, GL_UNSIGNED_BYTE, data);
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


// ***************** glGenBuffers *****************


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


// ***************** Other *****************


bool gen_mipmaps() {

	if (!glGenerateMipmap) return 0;
	glGenerateMipmap(GL_TEXTURE_2D);
	return 1;
}


