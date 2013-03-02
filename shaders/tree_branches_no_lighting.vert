varying vec3 normal;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position    = ftransform();
	normal         = gl_Normal; // world space (not normalized)
	gl_FrontColor  = gl_Color;
} 
