// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06

#ifndef _GL_EXT_ARB_H_
#define _GL_EXT_ARB_H_

#define MULTISAMPLE_ARB              0x809D
#define SAMPLE_ALPHA_TO_COVERAGE_ARB 0x809E
#define SAMPLE_ALPHA_TO_ONE_ARB      0x809F
#define SAMPLE_COVERAGE_ARB          0x80A0

#define MULTISAMPLE_BIT_ARB          0x20000000

#define SAMPLE_BUFFERS_ARB           0x80A8
#define SAMPLES_ARB                  0x80A9
#define SAMPLE_COVERAGE_VALUE_ARB    0x80AA
#define SAMPLE_COVERAGE_INVERT_ARB   0x80AB


void init_glew();


// multitexture prototypes
void setup_multitexture();
void set_multitex(unsigned tu_id);
void select_multitex(int id, unsigned tu_id, bool enable=1);
void disable_multitex(unsigned tu_id, bool reset);
void disable_multitex_a();
void multitex_coord2f(GLfloat s, GLfloat t, unsigned tu_id=0);
void multitex_coord2f_a(GLfloat s, GLfloat t);

// fog coord prototypes
void setup_fog_coord();
void set_fog_coord(GLfloat val);
void enable_fog_coord();
void disable_fog_coord();

// gen buffers prototypes
bool setup_gen_buffers();
unsigned create_vbo();
void bind_vbo(unsigned id, bool is_index=0);
void delete_vbo(unsigned id);
void upload_vbo_data(void const *const data, size_t size, bool is_index=0);
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index=0);
bool gen_mipmaps();


#endif // _GL_EXT_ARB_H_

