varying vec4 eye_space_pos;

void main()
{
	set_tc0_from_vert_id();
	gl_Position   = fg_ftransform();
	eye_space_pos = gl_ModelViewMatrix * fg_Vertex;
	gl_FrontColor = fg_Color;
}