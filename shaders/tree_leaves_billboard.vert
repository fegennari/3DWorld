varying vec4 eye_space_pos;
varying vec2 tc;

void main()
{
	tc            = gl_MultiTexCoord0;
	gl_Position   = ftransform();
	eye_space_pos = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor = gl_Color;
}