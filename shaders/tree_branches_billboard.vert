varying vec4 world_space_pos, eye_space_pos;

void main()
{
	set_tc0_from_vert_id();
	gl_Position     = ftransform();
	world_space_pos = gl_Vertex;
	eye_space_pos   = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor   = gl_Color;
} 
