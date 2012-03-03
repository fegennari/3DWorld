// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06

#ifndef _GL_EXT_ARB_H_
#define _GL_EXT_ARB_H_


void init_glew();


// multitexture prototypes
void setup_multitexture();
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
bool gen_mipmaps();


#endif // _GL_EXT_ARB_H_

