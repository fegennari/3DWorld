varying vec3 normal;
varying vec2 tc;

void main()
{
	tc            = gl_MultiTexCoord0;
	gl_Position   = ftransform();
	normal        = gl_Normal; // world space (not normalized)
	gl_FrontColor = gl_Color;
} 
