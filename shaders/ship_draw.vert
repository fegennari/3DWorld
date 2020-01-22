// for additional texture scaling control
uniform float tscale_s = 1.0;
uniform float tscale_t = 1.0;

out vec2 tc;
out vec4 epos;
out vec3 eye_norm;
#ifdef ENABLE_CUBE_MAP_REFLECT
uniform mat4 fg_ViewMatrix;
out vec3 vpos, ws_normal;
#endif

void main() {
	tc = fg_TexCoord*vec2(tscale_s, tscale_t);
	// Note: since we do a lot of transforms on the CPU and draw calls are small, it's more efficient to create the fg_NormalMatrix here
	eye_norm = normalize(transpose(inverse(mat3(fg_ModelViewMatrix))) * fg_Normal);
	epos     = fg_ModelViewMatrix * fg_Vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
	fg_Color_vf  = fg_Color;
	gl_PointSize = 1.0; // for particles, may be unnecessary
#ifdef ENABLE_CUBE_MAP_REFLECT
	ws_normal = normalize((transpose(fg_ViewMatrix) * vec4(eye_norm, 1)).xyz);
	vpos      = (inverse(fg_ViewMatrix) * epos).xyz; // world space
#endif
}
