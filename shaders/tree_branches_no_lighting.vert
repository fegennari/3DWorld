out vec3 normal;
out vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	gl_Position   = fg_ftransform();
	normal        = fg_Normal; // world space (not normalized)
	gl_FrontColor = fg_Color;
} 
