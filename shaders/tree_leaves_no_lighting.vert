varying vec3 normal;

void main()
{
	set_tc0_from_vert_id();
	gl_Position   = fg_ftransform();
	normal        = fg_Normal; // world space (not normalized)
	gl_FrontColor = fg_Color;
} 
