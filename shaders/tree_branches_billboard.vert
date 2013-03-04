varying vec4 world_space_pos, eye_space_pos;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0;
	gl_Position     = ftransform();
	world_space_pos = gl_Vertex;
	eye_space_pos   = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor   = gl_Color;
	gl_FogFragCoord = length(eye_space_pos.xyz);
} 
