// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06

#ifndef _GL_EXT_ARB_H_
#define _GL_EXT_ARB_H_


// multitexture prototypes
void setup_multitexture();
void set_multitex(unsigned tu_id);
void select_multitex(int id, unsigned tu_id);
void disable_multitex(unsigned tu_id, bool reset);
void disable_multitex_a();
void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id=0);
void multitex_coord2f_a(GLfloat s, GLfloat t);
bool has_multitex();

// fog coord prototypes
void setup_fog_coord_ext();
void set_fog_coord(GLfloat val);
void enable_fog_coord();
void disable_fog_coord();
bool has_fog_coord();

// gen buffers prototypes
bool setup_gen_buffers_arb();
unsigned create_vbo();
void bind_vbo(unsigned id, bool is_index=0);
void delete_vbo(unsigned id);
void upload_vbo_data(void const *const data, size_t size, bool is_index=0);
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index=0);


#endif // _GL_EXT_ARB_H_

