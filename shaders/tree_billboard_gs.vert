#ifdef USE_BINDLESS_TEXTURES
#extension GL_ARB_bindless_texture : require
layout (location = 4) in sampler2D color_handle; // 64-bit handle passed as two 32-bit uints
layout (location = 5) in sampler2D normal_handle;
#endif

out vec4 vertex_vs;
out vec4 color_vs;
out vec2 delta_vs;
#ifdef USE_BINDLESS_TEXTURES
out flat sampler2D color_tex_vs, normal_tex_vs;
#endif

void main() {
	vertex_vs = fg_Vertex;
	color_vs  = fg_Color;
	delta_vs  = fg_TexCoord.st; // Note: could use vec2 delta attribute
#ifdef USE_BINDLESS_TEXTURES
	color_tex_vs  = sampler2D(color_handle );
	normal_tex_vs = sampler2D(normal_handle);
#endif
}
