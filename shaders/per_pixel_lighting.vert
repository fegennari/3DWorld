uniform float point_size_pixels = 1.0; // for point sprite mode

varying vec4 epos;
varying vec3 normal;
varying vec2 tc;

void main()
{
	tc            = gl_MultiTexCoord0;
	normal        = normalize(gl_NormalMatrix * gl_Normal);
	epos          = gl_ModelViewMatrix * gl_Vertex;
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
#ifdef POINT_SPRITE_MODE
	gl_PointSize  = point_size_pixels;
#endif
}
