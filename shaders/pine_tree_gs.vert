void main()
{
	gl_FrontColor    = gl_Color;
	gl_Position      = gl_Vertex;
	gl_TexCoord[7].s = gl_Normal.x; // height comes from normal.x and is placed in TC7.s - FIXME
} 
