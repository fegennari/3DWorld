// for additional texture scaling control
uniform float tscale_s = 1.0;
uniform float tscale_t = 1.0;

out vec2 tc;
out vec4 epos;
out vec3 normal; // eye space

invariant gl_Position;

void main()
{
	tc = fg_TexCoord*vec2(tscale_s, tscale_t);
	// Note: since we do a lot of transforms on the CPU and draw calls are small, it's more efficient to create the fg_NormalMatrix here
	normal = normalize(transpose(inverse(mat3(fg_ModelViewMatrix))) * fg_Normal);
	epos   = fg_ModelViewMatrix * fg_Vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
	fg_Color_vf  = fg_Color;
	gl_PointSize = 1.0; // for particles, may be unnecessary
}
