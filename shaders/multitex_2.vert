out vec2 tc;

void main()
{
	tc          = fg_TexCoord;
	gl_Position = fg_ftransform();
	fg_Color_vf = fg_Color;
}
