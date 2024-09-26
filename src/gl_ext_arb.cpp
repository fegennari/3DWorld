// 3D World - GL EXT/ARB extension interface code
// by Frank Gennari
// 4/25/02
#include "3DWorld.h"
#include "gl_ext_arb.h"
#include "function_registry.h"
#include "inlines.h" // for render_to_texture_t::render() stuff


extern bool use_core_context;


void init_glew() {

// MacOSX check here, placeholder for eventual cross-platform porting
#if ((defined(__MACH__))&&(defined(__APPLE__)))
	// nothing
#else
	GLenum const err(glewInit());
	//glGetError(); // get and discard error flags (needed for unpatched GLEW)

	if (GLEW_OK != err) {
		std::cerr << "Error: " << glewGetErrorString(err) << endl;
		exit(1);
	}
#endif
	// Note: not sure why the 4.6 check is needed, but a user claimed that 4.6 was supported but 4.5 was not
	if (!glewIsSupported("GL_VERSION_4_5") && !glewIsSupported("GL_VERSION_4_6")) {
		std::cerr << "Error: GL version 4.5 not found; trying 4.4" << endl;

		if (!glewIsSupported("GL_VERSION_4_4")) { // does 3DWorld work with OpenGL 4.4? maybe most scenes?
			std::cerr << "Error: GL version 4.4 not found; exiting" << endl;
			exit(1);
		}
	}
}

GLenum get_internal_texture_format(int ncolors, bool compressed, bool linear_space) { // Note: ncolors=2 is unused
	assert(ncolors > 0 && ncolors <= 4);
	GLenum const cformats [4] = {GL_COMPRESSED_RED, GL_COMPRESSED_RG, GL_COMPRESSED_RGB_S3TC_DXT1_EXT, GL_COMPRESSED_RGBA_S3TC_DXT5_EXT}; // clouds and explosions look bad
	GLenum const formats  [4] = {GL_R8, GL_RG8, GL_RGB8, GL_RGBA8};
	// Note: only supports linear space for RGB and RGBA texture formats (others not checked)
	GLenum const scformats[4] = {GL_COMPRESSED_RED, GL_COMPRESSED_RG, GL_COMPRESSED_SRGB_S3TC_DXT1_EXT, GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT};
	GLenum const sformats [4] = {GL_R8, GL_RG8, GL_SRGB8, GL_SRGB8_ALPHA8};
	return (linear_space ? (compressed ? scformats : sformats) : (compressed ? cformats : formats))[ncolors-1];
}


void bind_texture_tu(unsigned tid, unsigned tu_id) {
	assert(tid);
	assert(tu_id < 32); // max is GL_TEXTURE31; too strict?
	glBindTextureUnit(tu_id, tid);
}


// ***************** 3D TEXTURES *****************


void bind_3d_texture(unsigned tid) {
	glBindTexture(GL_TEXTURE_3D, tid);
	assert(glIsTexture(tid));
}

void setup_3d_texture(unsigned &tid, int filter, int wrap) {

	glGenTextures(1, &tid);
	bind_3d_texture(tid);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, filter); // GL_LINEAR_MIPMAP_LINEAR?
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, filter);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, wrap);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, wrap);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, wrap);
}

unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data, int filter, int wrap, bool compress, unsigned bytes_per_pixel) {

	assert(data.size() == ncomp*bytes_per_pixel*xsz*ysz*zsz);
	unsigned tid(0);
	setup_3d_texture(tid, filter, wrap);

	if (bytes_per_pixel == 1) { // 8-bit texture
		glTexImage3D(GL_TEXTURE_3D, 0, get_internal_texture_format(ncomp, compress), xsz, ysz, zsz, 0, get_texture_format(ncomp), GL_UNSIGNED_BYTE, &data.front());
	}
	else if (bytes_per_pixel == 2) { // 16-bit texture
		assert(!compress);
		assert(ncomp > 0 && ncomp <= 4);
		GLenum const formats_16[4] = {GL_R16, GL_RG16, GL_RGB16, GL_RGBA16};
		glTexImage3D(GL_TEXTURE_3D, 0, formats_16[ncomp-1], xsz, ysz, zsz, 0, get_texture_format(ncomp), GL_UNSIGNED_SHORT, &data.front());
	}
	else {assert(0);} // not supported
	//gen_mipmaps(3);
	return tid;
}

void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data)
{
	bind_3d_texture(tid);
	glTexSubImage3D(GL_TEXTURE_3D, 0, xoff, yoff, zoff, xsz, ysz, zsz, get_texture_format(ncomp), GL_UNSIGNED_BYTE, data);
}


// ***************** VBOs / VAOs *****************


unsigned create_vbo() {
	unsigned vbo;
	assert(glGenBuffers);
	glGenBuffers(1, &vbo);
	assert(vbo > 0);
	return vbo;
}

int get_buffer_target(bool is_index) {return (is_index ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER);}

void bind_vbo(unsigned vbo, bool is_index) { // okay if vbo is zero
	if (use_core_context && vbo == 0) return; // no point in binding 0
	glBindBuffer(get_buffer_target(is_index), vbo);
	//if (vbo) {assert(glIsBuffer(vbo));}
}

void delete_vbo(unsigned vbo) {
	if (vbo != 0) {glDeleteBuffers(1, &vbo);}
}

unsigned mode_from_dynamic_level(int dynamic_level) {
	assert(dynamic_level <= 2);
	unsigned const modes[3] = {GL_STATIC_DRAW, GL_DYNAMIC_DRAW, GL_STREAM_DRAW};
	return modes[dynamic_level];
}
void upload_vbo_data(void const *const data, size_t size, bool is_index, int dynamic_level) {
	glBufferData(get_buffer_target(is_index), size, data, mode_from_dynamic_level(dynamic_level));
}

void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index) {
	glBufferSubData(get_buffer_target(is_index), offset, size, data);
}

void upload_vbo_sub_data_no_sync(void const *data, unsigned start_byte, unsigned size_bytes, bool is_index) {

	assert(data && size_bytes > 0);
	int const target(get_buffer_target(is_index));
	void *buffer(glMapBufferRange(target, start_byte, size_bytes, (GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_RANGE_BIT | GL_MAP_UNSYNCHRONIZED_BIT)));
	if (buffer == nullptr) {check_gl_error(132);}
	assert(buffer != nullptr);
	memcpy(buffer, data, size_bytes);
	glUnmapBuffer(target);
}

void bind_ubo(unsigned ubo) {
	glBindBuffer(GL_UNIFORM_BUFFER, ubo);
}
void bind_ubo_base(unsigned ubo, unsigned index) {
	glBindBufferBase(GL_UNIFORM_BUFFER, index, ubo);
}
void upload_ubo_sub_data(void const *const data, int offset, size_t size) {
	glBufferSubData(GL_UNIFORM_BUFFER, offset, size, data);
}
void upload_ubo_data(void const *const data, size_t size, int dynamic_level) {
	glBufferData(GL_UNIFORM_BUFFER, size, data, mode_from_dynamic_level(dynamic_level));
}

void ubo_wrap_t::allocate_with_size(unsigned size, int dynamic_level) {
	if (ubo_valid()) return; // already created, assumed at the correct size
	ubo = create_vbo();
	check_bind_ubo(ubo);
	upload_ubo_data(nullptr, size, dynamic_level);
	bind_ubo(0); // optional?
}


unsigned create_vao() {
	unsigned vao;
	assert(glGenVertexArrays);
	glGenVertexArrays(1, &vao); 
	assert(vao > 0);
	return vao;
}

unsigned default_vao(0);

void clear_default_vao() {delete_and_zero_vao(default_vao);}

void bind_vao(unsigned vao) { // okay if vao is zero
	if (use_core_context && vao == 0) {
		if (!default_vao) {default_vao = create_vao();}
		check_bind_vao(default_vao);
		return;
	}
	glBindVertexArray(vao);
}

void delete_vao(unsigned vao) {
	if (vao != 0) {glDeleteVertexArrays(1, &vao);}
}


void vbo_ring_buffer_t::ensure_vbo(unsigned min_size) {

	if (size < min_size) {
		clear_vbo(); // allocate a larger vbo
		size = max(2*size, min_size); // at least double
	}
	if (vbo) return; // done
	create_vbo_with_null_data(vbo, size, is_index, 1); // reserve the space but don't use it; use a dynamic draw buffer
	pos = 0;
}

void vbo_ring_buffer_t::orphan_buffer_and_reset() {
	bind_vbo(vbo, is_index);
	upload_vbo_data(NULL, size, is_index); // orphan the buffer (Note: fast on ATI, slow on Nvidia - pipeline flush?)
	bind_vbo(0, is_index); // not necessary?
	pos = 0; // wraparound to the beginning
}

void const *vbo_ring_buffer_t::add_verts_bind_vbo(void const *const v, unsigned size_bytes) {

	assert(v != NULL);
	assert(size_bytes > 0);
	ensure_vbo(4*size_bytes); // at least 4x the current data size
	if (!has_space_for(size_bytes)) {orphan_buffer_and_reset();} // end of buffer space
	bind_vbo(vbo, is_index);
	assert(has_space_for(size_bytes));
	upload_vbo_sub_data_no_sync(v, pos, size_bytes, is_index);
	void const *ret((unsigned char const *)((size_t)pos));
	pos += size_bytes; // data allocated
	align_vbo_ptr(pos); // 16-byte alignment - makes no difference?
	return ret;
}


// ***************** FBOs *****************

void bind_fbo_texture(unsigned fbo_id, unsigned tid, bool is_depth_fbo, bool multisample, bool is_array, unsigned *layer=nullptr) {
	assert(glIsTexture(tid));
	int const attachment(is_depth_fbo ? GL_DEPTH_ATTACHMENT : GL_COLOR_ATTACHMENT0);

	if (layer) {
		assert(!multisample); // untested; probably doesn't work
		glFramebufferTextureLayer(GL_FRAMEBUFFER, attachment, tid, 0, *layer);
	}
	else if (is_array) {
		assert(!multisample); // untested; probably doesn't work
		glFramebufferTexture(GL_FRAMEBUFFER, attachment, tid, 0);
	}
	else {
		glFramebufferTexture2D(GL_FRAMEBUFFER, attachment, get_2d_texture_target(0, multisample), tid, 0); // is_array=0
	}
	check_gl_error(551);
}

void create_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo, bool multisample, bool is_array, unsigned *layer) {
	// Create a framebuffer object
	check_gl_error(550);
	glGenFramebuffers(1, &fbo_id);
	bind_fbo(fbo_id);
	
	if (is_depth_fbo) { // Instruct openGL that we won't bind a color texture with the currently binded FBO
		glDrawBuffer(GL_NONE);
		glReadBuffer(GL_NONE);
	}
	// Attach the texture to FBO depth or color attachment point
	bind_fbo_texture(fbo_id, tid, is_depth_fbo, multisample, is_array, layer);
	// Check FBO status
	GLenum const status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
	assert(status != GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE);
	assert(status == GL_FRAMEBUFFER_COMPLETE);
}


void enable_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo, bool multisample, bool is_array, unsigned *layer) {
	if (!fbo_id) {create_fbo(fbo_id, tid, is_depth_fbo, multisample, is_array, layer);}
	assert(fbo_id > 0);
	bind_fbo(fbo_id); // rendering offscreen
}
void bind_fbo(unsigned fbo_id) {glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);}
void disable_fbo() {glBindFramebuffer(GL_FRAMEBUFFER, 0);}

void free_fbo(unsigned &fbo_id) {
	if (fbo_id > 0) {glDeleteFramebuffers(1, &fbo_id);}
	fbo_id = 0;
}

void bind_pbo(unsigned pbo_id) {glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id);}


unsigned create_depth_render_buffer(unsigned xsize, unsigned ysize, bool multisample) {

	unsigned depthrenderbuffer(0);
	glGenRenderbuffers(1, &depthrenderbuffer);
	assert(depthrenderbuffer > 0);
	glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
	if (multisample) {glRenderbufferStorageMultisample(GL_RENDERBUFFER, NUM_TEX_MS_SAMPLES, GL_DEPTH_COMPONENT16, xsize, ysize);}
	else {glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, xsize, ysize);}
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
	return depthrenderbuffer;
}
void disable_and_free_render_buffer(unsigned &render_buffer) {

	glBindRenderbuffer(GL_RENDERBUFFER, 0);
	if (render_buffer > 0) {glDeleteRenderbuffers(1, &render_buffer);}
	render_buffer = 0;
}


void render_to_texture_t::pre_render(float xsize, float ysize, unsigned nx, unsigned ny, point const &center, vector3d const &view_dir) const {

	assert(xsize > 0.0 && ysize > 0);
	assert(tsize > 0 && nx > 0 && ny > 0);
	// setup matrices
	glViewport(0, 0, nx*tsize, ny*tsize);
	fgMatrixMode(FG_PROJECTION);
	fgPushIdentityMatrix();
	fgOrtho(-xsize, xsize, -ysize, ysize, -(xsize + ysize), (xsize + ysize));
	fgMatrixMode(FG_MODELVIEW);
	fgPushIdentityMatrix();
	rotate_from_v2v(-plus_z, view_dir);
	translate_to(-center);
}


void render_to_texture_t::post_render() {

	restore_prev_mvm_pjm_state();
	//glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	disable_fbo();
	set_standard_viewport();
}


colorRGBA get_clear_color() {
	colorRGBA clear_color;
	glGetFloatv(GL_COLOR_CLEAR_VALUE, (float *)&clear_color.R);
	return clear_color;
}
void set_temp_clear_color(colorRGBA const &clear_color, bool clear_depth, bool clear_stencil) { // and also clear color, depth, and stencil buffers
	colorRGBA const orig_clear_color(get_clear_color());
	glClearColor_rgba(clear_color);
	glClear(GL_COLOR_BUFFER_BIT | (clear_depth ? GL_DEPTH_BUFFER_BIT : 0) | (clear_stencil ? GL_STENCIL_BUFFER_BIT : 0));
	glClearColor_rgba(orig_clear_color);
}


void render_to_texture_t::render(texture_pair_t &tpair, float xsize, float ysize, point const &center, vector3d const &view_dir,
	colorRGBA const &bkg_color, bool use_depth_buffer, bool mipmap)
{
	assert(!(mipmap && tpair.multisample));
	pre_render(xsize, ysize, 1, 1, center, view_dir); // setup matrices, etc.
	tpair.ensure_tid(tsize, mipmap);
	colorRGBA const clear_normal(0.5, 0.5, 0.5, 0.0);
	colorRGBA const clear_colors[2] = {bkg_color, clear_normal};
	unsigned fbo_id(0);
	enable_fbo(fbo_id, tpair.tids[0], 0, tpair.multisample); // too slow to create and free fbos every time?
	unsigned render_buffer(use_depth_buffer ? create_depth_render_buffer(tsize, tsize, tpair.multisample) : 0);

	for (unsigned d = 0; d < 2; ++d) { // {color, normal}
		if (d == 1) {bind_fbo_texture(fbo_id, tpair.tids[1], 0, tpair.multisample, 0);} // bind second texture; is_depth_fbo=0, is_array=0
		set_temp_clear_color(clear_colors[d], use_depth_buffer);
		draw_geom(d != 0);
		//if (tpair.multisample) {glBlitFramebuffer(...);} // ???
		if (mipmap) {build_texture_mipmaps(tpair.tids[d], 2);}
	}
	if (use_depth_buffer) {disable_and_free_render_buffer(render_buffer);}
	free_fbo(fbo_id);
	post_render(); // restore state
}


// ***************** SSBOs *****************

// See: https://www.geeks3d.com/20140704/tutorial-introduction-to-opengl-4-3-shader-storage-buffers-objects-ssbo-demo/
void bind_ssbo(unsigned ssbo) {glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);}

unsigned create_ssbo(unsigned data_sz, void const *const data) {
	unsigned ssbo;
	glGenBuffers(1, &ssbo);
	bind_ssbo(ssbo);
	if (data_sz > 0) {glBufferData(GL_SHADER_STORAGE_BUFFER, data_sz, data, GL_DYNAMIC_COPY);}
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	return ssbo;
}
void ensure_ssbo(unsigned &ssbo, unsigned data_sz, void const *const data) {
	if (ssbo == 0) {ssbo = create_ssbo(data_sz, data);}
}

void update_ssbo(unsigned ssbo, unsigned data_sz, void const *const data) {
	assert(ssbo != 0);
	assert(data_sz > 0 && data != nullptr);
	bind_ssbo(ssbo);
	void *p(glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY));
	assert(p != nullptr);
	memcpy(p, data, data_sz);
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
}

// unsigned block_index(glGetProgramResourceIndex(program, GL_SHADER_STORAGE_BLOCK, "shader_data"));
// glShaderStorageBlockBinding(program, block_index, 80);
// glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding_point_index, ssbo);


// ***************** Timers *****************

// see http://www.lighthouse3d.com/tutorials/opengl-short-tutorials/opengl-timer-query/
query_perf_timer_t::query_perf_timer_t () {glGenQueries   (1, &query_id);}
query_perf_timer_t::~query_perf_timer_t() {glDeleteQueries(1, &query_id);}

void query_perf_timer_t::time_query() { // records the time only after all previous commands have been completed
	glQueryCounter(query_id, GL_TIMESTAMP);
}
GLuint64 query_perf_timer_t::get_query() {
	int avail(0); // wait until the results are available
	while (!avail) {glGetQueryObjectiv(query_id, GL_QUERY_RESULT_AVAILABLE, &avail);}
	GLuint64 time;
	glGetQueryObjectui64v(query_id, GL_QUERY_RESULT, &time); // get query results
	return time;
}

gpu_timer_t::gpu_timer_t() {
	time_query.time_query();
	timer = get_timestamp();
}
void gpu_timer_t::show() {
	query_perf_timer_t time_query2;
	GLint64 timer2;
	time_query2.time_query();
	timer2 = get_timestamp();
	cout << "GPU time: " << 1.0E-6*(time_query2.get_query() - time_query.get_query()) << ", CPU time: " << 1.0E-6*(timer2 - timer) << endl;
}

GLint64 get_timestamp() {
	GLint64 timer;
	glGetInteger64v(GL_TIMESTAMP, &timer);
	return timer;
}


// ***************** Other *****************


bool gen_mipmaps(unsigned dim) { // cube maps = 6

	assert((dim >= 1 && dim <= 3) || dim == 6);
	if (!glGenerateMipmap) return 0;
	int const tex_dims[6] = {GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, 0, 0, GL_TEXTURE_CUBE_MAP};
	glGenerateMipmap(tex_dims[dim-1]);
	return 1;
}

void enable_instancing_for_shader_loc(int loc) {
	glEnableVertexAttribArray(loc);
	glVertexAttribDivisor(loc, 1);
}
void disable_instancing_for_shader_loc(int loc) {
	glVertexAttribDivisor(loc, 0);
	glDisableVertexAttribArray(loc);
}

bool check_for_tess_shader() {
	if (glPatchParameteri != nullptr) return 1; // found
	std::cerr << "*** Warning: Tessellation shader support not found; disabling tessellation" << endl;
	return 0;
}

void setup_stencil_buffer_write() {
	glClear(GL_STENCIL_BUFFER_BIT);
	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_ALWAYS, 0, ~0U);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // disable color writing, we only want to write to the stencil
	glDepthMask(GL_FALSE);
}
void end_stencil_write() {
	glDepthMask(GL_TRUE);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glDisable(GL_STENCIL_TEST);
}

GLuint DebugScope::global_scope_depth = 0;



