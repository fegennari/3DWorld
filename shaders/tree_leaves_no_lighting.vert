varying vec3 normal;

void main()
{
	set_tc0_from_vert_id();
	gl_Position   = ftransform();
	normal        = gl_Normal; // world space (not normalized)
	gl_FrontColor = gl_Color;
} 
