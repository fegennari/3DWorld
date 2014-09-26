out vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	gl_Position   = fg_ftransform();
	gl_FrontColor = fg_Color;
}
