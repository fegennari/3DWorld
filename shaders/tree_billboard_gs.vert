#ifdef USE_BINDLESS_TEXTURES
#extension GL_ARB_bindless_texture : require
layout (location = 4) in sampler2DArray tex_handle; // 64-bit handle passed as two 32-bit uints
#endif

out vec4 vertex_vs;
out vec4 color_vs;
out vec2 delta_vs;
#ifdef USE_BINDLESS_TEXTURES
out flat sampler2DArray color_normal_tex_vs;
#endif

void main() {
	vertex_vs = fg_Vertex;
	color_vs  = fg_Color;
	delta_vs  = fg_TexCoord.st; // Note: could use vec2 delta attribute
#ifdef USE_BINDLESS_TEXTURES
	color_normal_tex_vs = sampler2DArray(tex_handle);
#endif
}
