// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06

#ifndef _GL_EXT_ARB_H_
#define _GL_EXT_ARB_H_


inline GLenum get_internal_texture_format(int ncolors, bool compressed=0) {
	GLenum const cformats[4] = {GL_COMPRESSED_LUMINANCE, GL_COMPRESSED_LUMINANCE_ALPHA, GL_COMPRESSED_RGB, GL_COMPRESSED_RGBA};
	GLenum const formats [4] = {GL_LUMINANCE8, GL_LUMINANCE8_ALPHA8, GL_RGB8, GL_RGBA8};
	return (compressed ? cformats : formats)[ncolors-1];
}


inline GLenum get_texture_format(int ncolors) {
	assert(ncolors >= 1 && ncolors <= 4);
	GLenum const formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA}; // GL_BGRA is supposedly faster, but do we want to swap things here?
	return formats[ncolors-1];
}


void init_glew();


// multitexture prototypes
void setup_multitexture();
void set_active_texture(unsigned tu_id);
void set_multitex(unsigned tu_id);
void select_multitex(int id, unsigned tu_id, bool enable=1);
void disable_multitex(unsigned tu_id, bool reset);
void disable_multitex_a();
void multitex_coord_n(unsigned tu_id, float const *v, unsigned num);
void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id=0);
void multitex_coord2f_a(GLfloat s, GLfloat t);

// 3D texture prototypes
void bind_3d_texture(unsigned tid);
unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data, int filter, int wrap);
void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data);
void set_3d_texture_as_current(unsigned tid, unsigned tu_id);

// fog coord prototypes
void setup_fog_coord();
void set_fog_coord(GLfloat val);
void enable_fog_coord();
void disable_fog_coord();

// gl_ext_arb
bool setup_gen_buffers();
unsigned create_vbo();
void bind_vbo(unsigned id, bool is_index=0);
void delete_vbo(unsigned id);
void upload_vbo_data(void const *const data, size_t size, bool is_index=0);
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index=0);
void create_fbo(unsigned &fbo_id, unsigned depth_tid, bool is_depth_fbo);
void enable_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo);
void disable_fbo();
void free_fbo(unsigned &fbo_id);
unsigned create_depth_render_buffer(unsigned xsize, unsigned ysize);
void disable_and_free_render_buffer(unsigned &render_buffer);
bool gen_mipmaps(unsigned dim=2);


// templated vbo management utility functions/classes

template<typename T> void upload_to_vbo(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0) {
	assert(vbo > 0);
	bind_vbo(vbo, is_index);
	upload_vbo_data(&data.front(), data.size()*sizeof(T), is_index);
	if (end_with_bind0) {bind_vbo(0, is_index);}
}

template<typename T> bool create_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0) {
	if (vbo) return 0; // already uploaded
	vbo = create_vbo();
	upload_to_vbo(vbo, data, is_index, end_with_bind0);
	return 1;
}

template<typename T> void create_bind_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0) {
	if (!create_vbo_and_upload(vbo, data, is_index, 0)) {bind_vbo(vbo, is_index);}
}


struct indexed_vbo_manager_t {

	unsigned vbo, ivbo, gpu_mem;

	indexed_vbo_manager_t() : vbo(0), ivbo(0), gpu_mem(0) {}
	void reset_vbos_to_zero() {vbo = ivbo = gpu_mem = 0;};

	void clear_vbos() {
		delete_vbo(vbo);
		delete_vbo(ivbo);
		reset_vbos_to_zero();
	}
	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata) {
		if (!vbo ) {create_vbo_and_upload(vbo,  data,  0, 0); gpu_mem += data.size() *sizeof(vert_type_t );}
		if (!ivbo) {create_vbo_and_upload(ivbo, idata, 1, 0); gpu_mem += idata.size()*sizeof(index_type_t);}
	}
	void pre_render() const {
		assert(vbo && ivbo);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo, 1);
	}
	static void post_render() {
		bind_vbo(0, 0);
		bind_vbo(0, 1);
	}
};


#endif // _GL_EXT_ARB_H_

