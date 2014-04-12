uniform float point_size_pixels = 1.0; // for point sprite mode

varying vec4 epos;
varying vec3 normal;
varying vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	normal        = normalize(gl_NormalMatrix * fg_Normal);
	epos          = gl_ModelViewMatrix * fg_Vertex;
	gl_Position   = gl_ProjectionMatrix * epos;
	gl_FrontColor = fg_Color;
#ifdef POINT_SPRITE_MODE
	gl_PointSize  = point_size_pixels;
#endif
}
