// 3D World - GL EXT/ARB extension interface code
// by Frank Gennari
// 4/25/02
#include "GL/glew.h"
#include "3DWorld.h"
#define GL_GLEXT_PROTOTYPES 1
#include <GL/glext.h>
#include "gl_ext_arb.h"


void init_glew() {

	GLenum const err(glewInit());

	if (GLEW_OK != err) {
	  fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	  assert(0);
	}
	glEnable(GL_MULTISAMPLE_ARB);
}

void shader_test() {
	//glCreateShader(GL_VERTEX_SHADER);
	//glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
}


// ***************** MULTITEXTURING *****************

unsigned const MAX_MULTITEX = 32; // max is GL_TEXTURE31_ARB

unsigned max_used_multitex(0);
bool multitex_enabled[MAX_MULTITEX] = {0};

#define USE_MULTITEX // comment out to disable

#ifdef USE_MULTITEX

void set_multitex(unsigned tu_id) {

	static bool inited(0);
	if (!inited) setup_multitexture();
	inited = 1;
	assert(tu_id < MAX_MULTITEX); // Note: Assumes textures are defined sequentially
	multitex_enabled[tu_id] = 1;
	max_used_multitex = max(max_used_multitex, (tu_id+1));
	glActiveTextureARB(GL_TEXTURE0_ARB + tu_id);
}


void select_multitex(int id, unsigned tu_id) {

	set_multitex(tu_id);
	select_texture(id);
}


void disable_multitex(unsigned tu_id, bool reset) {

	assert(tu_id < max_used_multitex);
	multitex_enabled[tu_id] = 0;
	if ((tu_id+1) == max_used_multitex) --max_used_multitex;
	glActiveTextureARB((GL_TEXTURE0_ARB + tu_id));
	glDisable(GL_TEXTURE_2D);
	if (reset) glActiveTextureARB(GL_TEXTURE0_ARB); // end back at texture 0
}


void disable_multitex_a() {

	for (unsigned i = 0; i < max_used_multitex; ++i) {
		if (multitex_enabled[i]) disable_multitex(i, 0);
	}
	glActiveTextureARB(GL_TEXTURE0_ARB); // end back at texture 0
	max_used_multitex = 0;
}


void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id) {

	assert(tu_id < max_used_multitex);
	glMultiTexCoord2fARB((GL_TEXTURE0_ARB + tu_id), s, t);
	//glTexCoord2f(s, t);
}


void multitex_coord2f_a(GLfloat s, GLfloat t) {

	for (unsigned i = 0; i < max_used_multitex; ++i) {
		if (multitex_enabled[i]) multitex_coord2f(s, t, i);
	}
}


void setup_multitexture() { // Windows specific

	if (!glMultiTexCoord2fARB) {
		cout << "*** Can't find GL_ARB_multitexture extension ***" << endl;
		assert(0);
	}
	//glClientActiveTextureARB(GL_TEXTURE0_ARB);
	//glTexCoordPointer(2, GL_FLOAT, 0, tp0);
	//glEnableClientState(GL_TEXTURE_COORD_ARRAY);
}


bool has_multitex() {return 1;}

#else

void set_multitex(unsigned tu_id) {}
void select_multitex(int id, unsigned tu_id) {}
void disable_multitex(unsigned tu_id, bool reset) {}
void disable_multitex_a() {}
void multitex_coord2f(float s, float t, unsigned tu_id) {}
void multitex_coord2f_a(float s, float t) {}
void setup_multitexture() {}
bool has_multitex() {return 0;}

#endif // USE_MULTITEX


// ***************** FOG_COORD_EXT *****************

#define USE_FOG_COORD_EXT // comment out to disable

#ifdef USE_FOG_COORD_EXT


void setup_fog_coord_ext() {

	if (!glFogCoordfEXT) {
		cout << "*** Can't find EXT_fog_coord extension ***" << endl;
		assert(0);
	}
}


void set_fog_coord(GLfloat val) {
	glFogCoordfEXT(val);
}

void enable_fog_coord() {
	glFogi(GL_FOG_COORDINATE_SOURCE_EXT, GL_FOG_COORDINATE_EXT);
}

void disable_fog_coord() {
	glFogi(GL_FOG_COORDINATE_SOURCE_EXT, GL_FRAGMENT_DEPTH_EXT); // is this correct?
}

bool has_fog_coord() {return (glFogCoordfEXT != NULL);}

#else

void setup_fog_coord_ext() {}
void set_fog_coord(float val) {}
void enable_fog_coord() {}
void disable_fog_coord() {}
bool has_fog_coord() {return 0;}

#endif // USE_FOG_COORD_EXT


// ***************** glGenBuffersARB *****************

#define USE_VBOS // comment out to disable

#ifdef USE_VBOS


bool setup_gen_buffers_arb() {

	static int retval(2);
	
	if (retval == 2) { // not yet setup
		if (!glGenBuffersARB) {
			cout << "*** glGenBuffersARB is not supported ***" << endl;
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
	glGenBuffersARB(1, &id);
	return id;
}

int get_buffer_target(bool is_index) {return (is_index ? GL_ELEMENT_ARRAY_BUFFER_ARB : GL_ARRAY_BUFFER_ARB);}

void bind_vbo(unsigned id, bool is_index) {
	glBindBufferARB(get_buffer_target(is_index), id);
}

void delete_vbo(unsigned id) {
	if (id == 0) return;
	glDeleteBuffersARB(1, &id);
}

void upload_vbo_data(void const *const data, size_t size, bool is_index) {
	// hard coded for drawing of static data
	glBufferDataARB(get_buffer_target(is_index), size, data, GL_STATIC_DRAW_ARB);
}

void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index) {
	glBufferSubDataARB(get_buffer_target(is_index), offset, size, data);
}

#else

bool setup_gen_buffers_arb() {return 0;}
unsigned create_vbo() {return 0;}
void bind_vbo(unsigned id, bool is_index) {}
void delete_vbo(unsigned id) {}
void upload_vbo_data(void const *const data, size_t size, bool is_index) {}
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index) {}

#endif // USE_VBOS


