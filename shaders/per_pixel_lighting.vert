uniform float point_size_pixels   = 1.0; // for point sprite mode
uniform float vertex_offset_scale = 0.0; // hack to make vertex_offset ignored when unused/unset

in vec3 vertex_offset; // not always used

out vec4 epos;
out vec3 normal; // eye space
out vec2 tc;

void main() {
	tc          = fg_TexCoord;
	normal      = normalize(fg_NormalMatrix * fg_Normal);
	vec4 vertex = vec4(vertex_offset_scale*vertex_offset, 0.0) + fg_Vertex;
	epos        = fg_ModelViewMatrix * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	fg_Color_vf = fg_Color;
#ifdef POINT_SPRITE_MODE
	gl_PointSize = point_size_pixels;
#endif
}
