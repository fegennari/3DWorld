uniform float point_size_pixels = 1.0; // for point sprite mode

out vec4 epos;
out vec3 normal;
out vec2 tc;

void main()
{
	tc          = fg_TexCoord;
	normal      = normalize(fg_NormalMatrix * fg_Normal);
	epos        = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	fg_Color_vf = fg_Color;
#ifdef POINT_SPRITE_MODE
	gl_PointSize = point_size_pixels;
#endif
}
