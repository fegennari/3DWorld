varying vec4 world_space_pos, eye_space_pos;
varying vec2 tc;

void main()
{
	tc              = gl_MultiTexCoord0;
	gl_Position     = ftransform();
	world_space_pos = gl_Vertex;
	eye_space_pos   = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor   = gl_Color;
} 
