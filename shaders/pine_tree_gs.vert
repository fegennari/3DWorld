void main()
{
	gl_FrontColor    = fg_Color;
	gl_Position      = fg_Vertex;
	gl_TexCoord[7].s = fg_Normal.x; // height comes from normal.x and is placed in TC7.s - FIXME
} 
