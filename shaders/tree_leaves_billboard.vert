varying vec4 eye_space_pos;

void main()
{
	set_tc0_from_vert_id();
	gl_Position   = ftransform();
	eye_space_pos = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor = gl_Color;
}